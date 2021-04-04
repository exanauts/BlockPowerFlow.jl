module CUSOLVERRF

using LinearAlgebra
using SparseArrays
using CUDA
using CUDA.CUSPARSE
using CUDA.CUSOLVER

import CUDA.CUBLAS: unsafe_batch

include("libcusolverRf_common.jl")
include("libcusolverRf_api.jl")

function cusolverRfCreate()
    handle_ref = Ref{cusolverRfHandle_t}()
    cusolverRfCreate(handle_ref)
    return handle_ref[]
end

function _cusolverRfDestroy(handle)
    if handle != C_NULL
        cusolverRfDestroy(handle)
        handle = C_NULL
    end
end

abstract type CusolverFactorization{T} <: LinearAlgebra.Factorization{T} end

struct CusolverRfLU{T} <: CusolverFactorization{T}
    gH::Ptr{cusolverRfHandle_t}
    nrhs::Int
    n::Int
    m::Int
    nnzA::Int
    drowsA::CuVector{Cint}
    dcolsA::CuVector{Cint}
    dP::CuVector{Cint}
    dQ::CuVector{Cint}
    dT::CuVector{T}
end

# Batch mode should not mix with classical LU factorization.
# We implement a structure apart.
struct CusolverRfLUBatch{T} <: CusolverFactorization{T}
    gH::Ptr{cusolverRfHandle_t}
    batchsize::Int
    n::Int
    m::Int
    nnzA::Int
    drowsA::CuVector{Cint}
    dcolsA::CuVector{Cint}
    dP::CuVector{Cint}
    dQ::CuVector{Cint}
    dT::CuVector{T}
end

cudestroy!(rflu::CusolverFactorization) = _cusolverRfDestroy(rflu.gH)

function glu_analysis_host(
    n, m, nnzA, h_rowsA, h_colsA, h_valsA,
    ordering, tol,
)
    h_Qreorder = zeros(Cint, n)

    spH = CUSOLVER.sparse_handle()

    # Create matrix descriptor
    desca = CUSPARSE.CuMatrixDescriptor()
    CUSPARSE.cusparseSetMatType(desca, CUSPARSE.CUSPARSE_MATRIX_TYPE_GENERAL)
    CUSPARSE.cusparseSetMatIndexBase(desca, CUSPARSE.CUSPARSE_INDEX_BASE_ZERO)

    # Reordering
    if ordering == :AMD
        CUSOLVER.cusolverSpXcsrsymamdHost(
            spH,
            n, nnzA, desca,
            h_rowsA, h_colsA, h_Qreorder,
        )
    elseif ordering == :MDQ
        CUSOLVER.cusolverSpXcsrsymmdqHost(
            spH,
            n, nnzA, desca,
            h_rowsA, h_colsA, h_Qreorder,
        )

    elseif ordering == :METIS
        CUSOLVER.cusolverSpXcsrmetisndHost(
            spH,
            n, nnzA, desca,
            h_rowsA, h_colsA, C_NULL, h_Qreorder,
        )
    elseif ordering == :RCM
        CUSOLVER.cusolverSpXcsrsymrcmHost(
            spH,
            n, nnzA, desca,
            h_rowsA, h_colsA, h_Qreorder,
        )
    end

    # Create duplicate matrix for reordering
    h_rowsB = copy(h_rowsA)
    h_colsB = copy(h_colsA)
    h_valsB = copy(h_valsA)
    h_mapBfromA = zeros(Cint, nnzA)
    @inbounds for i in 1:nnzA
        h_mapBfromA[i] = i # identity matrix
    end

    # Compute permutation in two steps
    size_perm = Ref{Csize_t}(0)
    CUSOLVER.cusolverSpXcsrperm_bufferSizeHost(
        spH,
        m, n, nnzA, desca,
        h_rowsB, h_colsB, h_Qreorder, h_Qreorder,
        size_perm,
    )

    buffer_cpu = Base.Libc.malloc(size_perm[] * sizeof(Cint))
    status = CUSOLVER.cusolverSpXcsrpermHost(
        spH,
        m, n, nnzA, desca,
        h_rowsB, h_colsB, h_Qreorder, h_Qreorder, h_mapBfromA,
        buffer_cpu,
    )

    # Apply permutation
    h_valsB = h_valsA[h_mapBfromA]

    # LU Factorization
    info = Ref{CUSOLVER.csrqrInfo_t}()
    status = cusolverSpCreateCsrluInfoHost(info)

    status = cusolverSpXcsrluAnalysisHost(
        spH,
        m, nnzA, desca,
        h_rowsB, h_colsB, info[],
    )

    size_internal = Ref{Cint}(0)
    size_lu = Ref{Cint}(0)
    status = cusolverSpDcsrluBufferInfoHost(
        spH,
        n, nnzA, desca,
        h_valsB, h_rowsB, h_colsB,
        info[],
        size_internal, size_lu
    )

    n_bytes = size_lu[] * sizeof(Cint)
    println("Run LU factorization on the host with ", n_bytes, " bytes")
    buffer_lu = Base.Libc.malloc(n_bytes)
    pivot_threshold = 1.0

    status = cusolverSpDcsrluFactorHost(
        spH, n, nnzA, desca,
        h_valsB, h_rowsB, h_colsB,
        info[], pivot_threshold,
        buffer_lu,
    )

    # Check singularity
    singularity = Ref{Cint}(0)
    status = cusolverSpDcsrluZeroPivotHost(
        spH, info[], tol, singularity,
    )

    # Check that the matrix is nonsingular
    if singularity[] >= 0
        SingularException(singularity[])
    end

    # Get size of L and U
    pnnzU = Ref{Cint}(0)
    pnnzL = Ref{Cint}(0)
    status = cusolverSpXcsrluNnzHost(
        spH,
        pnnzL, pnnzU, info[],
    )

    nnzL = pnnzL[]
    nnzU = pnnzU[]

    # Retrieve L and U matrices
    h_Plu = zeros(Cint, m)
    h_Qlu = zeros(Cint, n)

    h_valsL = zeros(nnzL)
    h_rowsL = zeros(Cint, m+1)
    h_colsL = zeros(Cint, nnzL)

    h_valsU = zeros(nnzU)
    h_rowsU = zeros(Cint, m+1)
    h_colsU = zeros(Cint, nnzU)

    # Extract
    status = cusolverSpDcsrluExtractHost(
        spH,
        h_Plu, h_Qlu,
        desca,
        h_valsL, h_rowsL, h_colsL,
        desca,
        h_valsU, h_rowsU, h_colsU,
        info[],
        buffer_lu,
    )

    h_P = h_Qreorder[h_Plu .+ 1]
    h_Q = h_Qreorder[h_Qlu .+ 1]

    if buffer_lu != C_NULL
        Base.Libc.free(buffer_lu)
        buffer_lu = C_NULL
    end
    if buffer_cpu != C_NULL
        Base.Libc.free(buffer_cpu)
        buffer_cpu = C_NULL
    end

    return (
        nnzL=nnzL, rowsL=h_rowsL, colsL=h_colsL, valsL=h_valsL,
        nnzU=nnzU, rowsU=h_rowsU, colsU=h_colsU, valsU=h_valsU,
        P=h_P, Q=h_Q,
    )
end

function CusolverRfLU(
    A::AbstractCuSparseMatrix{T};
    nrhs=1, ordering=:AMD, tol=1e-8, fast_mode=true,
) where T
    m, n = size(A)
    @assert m == n # only squared matrices are supported
    if nrhs > 1
        error("Currently CusolverRF supports only one right-hand side.")
    end

    nnzA = nnz(A)
    # Allocations (device)
    d_P = CUDA.zeros(Cint, m)
    d_Q = CUDA.zeros(Cint, n)
    d_T = CUDA.zeros(Cdouble, m * nrhs)

    # Allocations (host)
    # Transfer data to host
    if isa(A, CuSparseMatrixCSR)
        h_rowsA = A.rowPtr |> Vector{Cint}
        h_colsA = A.colVal |> Vector{Cint}
        h_valsA = A.nzVal |> Vector{T}
    elseif isa(A, CuSparseMatrixCSC)
        h_rowsA = A.colPtr |> Vector{Cint}
        h_colsA = A.rowVal |> Vector{Cint}
        h_valsA = A.nzVal |> Vector{T}
    end

    # cusolverRf is 0-based
    h_rowsA .-= 1
    h_colsA .-= 1

    lu = glu_analysis_host(n, m, nnzA, h_rowsA, h_colsA, h_valsA, ordering, tol)

    # Create handle
    gH = cusolverRfCreate()

    ## OPTIONS
    # Set fast mode
    if fast_mode
        status = cusolverRfSetResetValuesFastMode(gH, CUSOLVERRF_RESET_VALUES_FAST_MODE_ON)
    else
        status = cusolverRfSetResetValuesFastMode(gH, CUSOLVERRF_RESET_VALUES_FAST_MODE_OFF)
    end
    # Set numeric properties
    nzero = 0.0
    nboost = 0.0
    status = cusolverRfSetNumericProperties(gH, nzero, nboost)

    # Set matrix format
    status = cusolverRfSetMatrixFormat(
        gH,
        CUSOLVERRF_MATRIX_FORMAT_CSR,
        CUSOLVERRF_UNIT_DIAGONAL_ASSUMED_L
    )

    status = cusolverRfSetAlgs(
        gH,
        CUSOLVERRF_FACTORIZATION_ALG2,
        CUSOLVERRF_TRIANGULAR_SOLVE_ALG2,
    )

    # Assemble internal data structures
    status = cusolverRfSetupHost(
        n, nnzA, h_rowsA, h_colsA, h_valsA,
        lu.nnzL, lu.rowsL, lu.colsL, lu.valsL,
        lu.nnzU, lu.rowsU, lu.colsU, lu.valsU,
        lu.P, lu.Q,
        gH
    )
    # Analyze available parallelism
    status = CUDA.@sync cusolverRfAnalyze(gH)
    # LU refactorization
    status = CUDA.@sync cusolverRfRefactor(gH)

    return CusolverRfLU{T}(
        gH, nrhs, n, m, nnzA,
        h_rowsA, h_colsA, lu.P, lu.Q, d_T
    )
end

# Update factorization inplace
function update!(rflu::CusolverRfLU{T}, A::AbstractCuSparseMatrix{T}) where T
    status = CUDA.@sync cusolverRfResetValues(
        rflu.n, rflu.nnzA,
        rflu.drowsA, rflu.dcolsA, A.nzVal, rflu.dP, rflu.dQ,
        rflu.gH
    )

    # LU refactorization
    status = CUDA.@sync cusolverRfRefactor(rflu.gH)
    return
end

# Solve linear-system
function solve!(x::CuArray{T}, rflu::CusolverRfLU{T}, b::CuArray{T}) where T
    n = rflu.n
    nrhs = 1
    copyto!(x, b)
    # Forward and backward solve
    status = CUDA.@sync cusolverRfSolve(rflu.gH, rflu.dP, rflu.dQ, nrhs, rflu.dT, n, x, n)
    return
end

LinearAlgebra.lu(A::AbstractCuSparseMatrix; options...) = CusolverRfLU(A; options...)
LinearAlgebra.lu!(rflu::CusolverRfLU, A::AbstractCuSparseMatrix) = update!(rflu, A)
LinearAlgebra.ldiv!(x::CuArray, rflu::CusolverRfLU, b::CuArray) = solve!(x, rflu, b)

function CusolverRfLUBatch(
    A::CuSparseMatrixCSR{T}, batchsize::Int;
    ordering=:AMD, tol=1e-8, fast_mode=false,
) where T
    m, n = size(A)
    @assert m == n # only squared matrices are supported

    if batchsize == 1
        error("Please specify a number of batch greater than 1 for `CuSolverRfLUBatch`")
    end
    nnzA = nnz(A)
    # Allocations (device)
    d_P = CUDA.zeros(Cint, m)
    d_Q = CUDA.zeros(Cint, n)
    d_T = CUDA.zeros(Cdouble, m * batchsize * 2)

    # Allocations (host)
    # Transfer data to host
    h_rowsA = A.rowPtr |> Vector{Cint}
    h_colsA = A.colVal |> Vector{Cint}
    h_valsA = A.nzVal |> Vector{T}

    # cusolverRf is 0-based
    h_rowsA .-= 1
    h_colsA .-= 1

    lu = glu_analysis_host(n, m, nnzA, h_rowsA, h_colsA, h_valsA, ordering, tol)

    # Create handle
    gH = cusolverRfCreate()

    ## OPTIONS
    # Set fast mode
    if fast_mode
        status = cusolverRfSetResetValuesFastMode(gH, CUSOLVERRF_RESET_VALUES_FAST_MODE_ON)
    else
        status = cusolverRfSetResetValuesFastMode(gH, CUSOLVERRF_RESET_VALUES_FAST_MODE_OFF)
    end
    # Set numeric properties
    nzero = 0.0
    nboost = 0.0
    status = cusolverRfSetNumericProperties(gH, nzero, nboost)

    # Set matrix format
    status = cusolverRfSetMatrixFormat(
        gH,
        CUSOLVERRF_MATRIX_FORMAT_CSR,
        CUSOLVERRF_UNIT_DIAGONAL_ASSUMED_L
    )

    status = cusolverRfSetAlgs(
        gH,
        CUSOLVERRF_FACTORIZATION_ALG2,
        CUSOLVERRF_TRIANGULAR_SOLVE_ALG2,
    )

    # Assemble internal data structures
    h_valsA_batch = zeros(batchsize, nnzA)
    h_valsA_batch = Vector{Float64}[copy(h_valsA) for i in 1:batchsize]
    ptrA_batch = pointer.(h_valsA_batch)
    status = cusolverRfBatchSetupHost(
        batchsize,
        n, nnzA, h_rowsA, h_colsA, ptrA_batch,
        lu.nnzL, lu.rowsL, lu.colsL, lu.valsL,
        lu.nnzU, lu.rowsU, lu.colsU, lu.valsU,
        lu.P, lu.Q,
        gH
    )
    # Analyze available parallelism
    status = CUDA.@sync cusolverRfBatchAnalyze(gH)
    # LU refactorization
    status = CUDA.@sync cusolverRfBatchRefactor(gH)

    return CusolverRfLUBatch{T}(
        gH, batchsize, n, m, nnzA,
        h_rowsA, h_colsA, lu.P, lu.Q, d_T
    )
end

# Update factorization inplace
function update!(rflu::CusolverRfLUBatch{T}, A::CuSparseMatrixCSR{T}) where T
    status = CUDA.@sync cusolverRfResetValues(
        rflu.n, rflu.nnzA,
        rflu.drowsA, rflu.dcolsA, A.nzVal, rflu.dP, rflu.dQ,
        rflu.gH
    )

    # LU refactorization
    status = CUDA.@sync cusolverRfRefactor(rflu.gH)
    return
end

# Solve linear-system
function solve!(x::Vector{CuVector{T}}, rflu::CusolverRfLUBatch{T}, b::Vector{CuVector{T}}) where T
    @assert rflu.batchsize > 1
    n = rflu.n
    nrhs = 1
    for i in 1:rflu.batchsize
        copyto!(x[i], b[i])
    end
    Xptrs = unsafe_batch(x)
    # Forward and backward solve
    status = CUDA.@sync cusolverRfBatchSolve(rflu.gH, rflu.dP, rflu.dQ, nrhs, rflu.dT, n, Xptrs, n)
    return
end

end

