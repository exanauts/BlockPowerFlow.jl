module CUSOLVERRF

using SparseArrays
using CUDA
using CUDA.CUSPARSE
using CUDA.CUSOLVER

include("libcusolverRf_common.jl")
include("libcusolverRf_api.jl")

function cusolverRfCreate()
  handle_ref = Ref{cusolverRfHandle_t}()
  cusolverRfCreate(handle_ref)
  return handle_ref[]
end

function _cusolverRfDestroy(handle)
  handle_ref = Ref{cusolverRfHandle_t}(handle)
  cusolverRfCreate(handle_ref)
end

function batch_lu!(A::CUSPARSE.CuSparseMatrixCSR{Float64},
                   b::CUDA.CuVector{Float64},
                   x::CUDA.CuVector{Float64},
                   tol::Float64,
                   batchsize::Cint,
                   inda::Char)
    m, n = size(A)
    @assert m == n

    nnzA = nnz(A)
    nrhs = 1 # currently only 1 is supported

    # Allocations (host)
    h_P = zeros(m)
    h_Q = zeros(n)

    # Allocations (device)
    d_P = CUDA.zeros(Cint, m)
    d_Q = CUDA.zeros(Cint, n)
    d_T = CUDA.zeros(Cdouble, m * nrhs)

    # Load data
    d_rowsA = A.rowPtr
    d_colsA = A.colVal
    d_valsA = A.nzVal

    # Transfer data to host
    h_rowsA = d_rowsA |> Vector{Cint}
    h_colsA = d_colsA |> Vector{Cint}
    h_valsA = d_valsA |> Vector{Float64}

    # cusolverRf is base-0
    h_rowsA .-= 1
    h_colsA .-= 1

    h_Qreorder = zeros(Cint, n)

    spH = CUSOLVER.sparse_handle()
    # First factorization using GLU
    desca = CUSPARSE.CuMatrixDescriptor()
    CUSPARSE.cusparseSetMatType(desca, CUSPARSE.CUSPARSE_MATRIX_TYPE_GENERAL)
    CUSPARSE.cusparseSetMatIndexBase(desca, CUSPARSE.CUSPARSE_INDEX_BASE_ZERO)

    # reordering
    CUSOLVER.cusolverSpXcsrsymrcmHost(
        spH,
        n, nnzA, desca,
        h_rowsA, h_colsA, h_Qreorder,
    )
    h_rowsB = copy(h_rowsA)
    h_colsB = copy(h_colsA)
    h_valsB = copy(h_valsA)

    size_perm = Ref{Csize_t}(0)
    CUSOLVER.cusolverSpXcsrperm_bufferSizeHost(
        spH,
        m, n, nnzA, desca,
        h_rowsB, h_colsB, h_Qreorder, h_Qreorder,
        size_perm,
    )

    h_mapBfromA = zeros(Cint, nnzA)
    for i in 1:nnzA
        h_mapBfromA[i] = i
    end

    buffer_cpu = Base.Libc.malloc(size_perm[] * sizeof(Cint))
    status = CUSOLVER.cusolverSpXcsrpermHost(
        spH,
        m, n, nnzA, desca,
        h_rowsB, h_colsB, h_Qreorder, h_Qreorder, h_mapBfromA,
        buffer_cpu,
    )

    # Apply permutation
    h_valsB = h_valsA[h_mapBfromA]

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

    CUDA.@time status = cusolverSpDcsrluFactorHost(
        spH, n, nnzA, desca,
        h_valsB, h_rowsB, h_colsB,
        info[], pivot_threshold,
        buffer_lu,
    )

    # Check singularity
    singularity = Ref{Cint}(0)
    status = CUDA.@time cusolverSpDcsrluZeroPivotHost(
        spH, info[], tol, singularity,
    )

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
    status = CUDA.@time cusolverSpDcsrluExtractHost(
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

    ##################################
    # Create handle
    gH = cusolverRfCreate()

    ## OPTIONS
    # Set fast mode
    status = cusolverRfSetResetValuesFastMode(gH, CUSOLVERRF_RESET_VALUES_FAST_MODE_ON)
    # Set options
    nzero = 0.0
    nboost = 0.0
    status = cusolverRfSetNumericProperties(gH, nzero, nboost)
    # Set matrix format
    status = cusolverRfSetMatrixFormat(gH, CUSOLVERRF_MATRIX_FORMAT_CSR, CUSOLVERRF_UNIT_DIAGONAL_ASSUMED_L)
    status = cusolverRfSetAlgs(
        gH,
        CUSOLVERRF_FACTORIZATION_ALG2,
        CUSOLVERRF_TRIANGULAR_SOLVE_ALG2,
    )

    # Assemble internal data structures
    status = cusolverRfSetupHost(
        n, nnzA, h_rowsA, h_colsA, h_valsA,
        nnzL, h_rowsL, h_colsL, h_valsL,
        nnzU, h_rowsU, h_colsU, h_valsU,
        h_P, h_Q,
        gH
    )
    # TODO: add cudaSynchronise?

    # Analyze available parallelism
    status = CUDA.@time cusolverRfAnalyze(gH)

    d_rowsA = CUDA.zeros(Cint, n+1)
    d_colsA = CUDA.zeros(Cint, nnzA)
    d_valsA = CUDA.zeros(Cdouble, nnzA)
    copy!(d_rowsA, h_rowsA)
    copy!(d_colsA, h_colsA)
    copy!(d_valsA, h_valsA)
    copy!(d_P, h_P)
    copy!(d_Q, h_Q)

    status = CUDA.@time cusolverRfResetValues(
        n, nnzA,
        d_rowsA, d_colsA, d_valsA, d_P, d_Q,
        gH
    )

    # LU re-factorization
    status = CUDA.@time cusolverRfRefactor(gH)

    copy!(x, b)
    # Forward and backward solve
    status = CUDA.@time cusolverRfSolve(gH, d_P, d_Q, nrhs, d_T, n, x, n)
    status = @time cusolverRfSolve(gH, d_P, d_Q, nrhs, d_T, n, x, n)
    copy!(x, b)
    status = @time cusolverRfSolve(gH, d_P, d_Q, nrhs, d_T, n, x, n)

    _cusolverRfDestroy(gH)
    println(gH)

    return x
end

end
