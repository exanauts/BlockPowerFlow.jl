using MatrixMarket
using SparseArrays
using LinearAlgebra

import KernelAbstractions
using CUDA
using CUDA.CUSPARSE

import ExaPF
import BlockPowerFlow

const LS = ExaPF.LinearSolvers
const BPF = BlockPowerFlow

# name = "case9"
# name = "case57"
# name = "case118"
# name = "case300"
# name = "case1354pegase"
# name = "case_ACTIVSg2000"
# name = "case2869pegase"
name = "case9241pegase"
# name = "case13659pegase"
# name = "case_ACTIVSg10k"
# name = "case_ACTIVSg25k"
# name = "case_ACTIVSg70k"

if CUDA.functional()
    path_A = joinpath(@__DIR__, "dump", "Gx_$(name).mtx")
    path_B = joinpath(@__DIR__, "dump", "Gu_$(name).mtx")

    # Storage on CPU
    A = mmread(path_A)
    B = mmread(path_B)
    B = Matrix(B)

    # Storage on GPU
    A = CuSparseMatrixCSR(A)
    B = CuMatrix(B)

    # Create the ILU(0) preconditioner
    m,n = size(A)
    s,p = size(B)
    P = BPF.ilu0(A, p)

    # Create a block-Jacobi preconditioner
    size_block = 32
    nblocks = div(n, size_block)
    device = ExaPF.CUDABackend()
    BJ = LS.BlockJacobiPreconditioner(A, nblocks, device)
    LS.update(BJ, A, device)

    # Parameters
    verbose = 1
    rtol = 0.0
    atol = 1e-10
    memory = 1
    restart = false
    itmax = 50

    println("Problem $name of size ($m,$n) with $p right-hand sides.\n")

    # Solve the linear system with a right preconditioner ILU(0)
    println("Problem $name with a right preconditioner ILU(0)")
    X_gmres, stats = BPF.block_gmres(A, B; N=P, ldiv=true, verbose, atol, rtol, memory, restart, itmax)
    RNorm = norm(B - A * X_gmres)
    println("‖Rₖ‖: ", RNorm)

    # Solve the linear system with a right preconditioner block-Jacobi
    println("Problem $name with a right preconditioner block-Jacobi")
    X_gmres, stats = BPF.block_gmres(A, B; N=BJ, verbose, atol, rtol, memory, restart, itmax)
    RNorm = norm(B - A * X_gmres)
    println("‖Rₖ‖: ", RNorm)

    # Solve the linear system without a preconditioner
    println("Problem $name without a preconditioner")
    X_gmres, stats = BPF.block_gmres(A, B; verbose, atol, rtol, memory, restart, itmax)
    RNorm = norm(B - A * X_gmres)
    println("‖Rₖ‖: ", RNorm)
end
