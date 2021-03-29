using Test
using LinearAlgebra
using MatrixMarket

using BlockPowerFlow

using CUDA
using CUDA.CUSPARSE
import ExaPF
import ExaPF: LinearSolvers

const LS = LinearSolvers
const BPF = BlockPowerFlow

function blk_bicgstab(J; nrhs=10, nblocks=8, gpu=false)
    n = size(J, 1)
    # Init preconditioner
    precond = LS.BlockJacobiPreconditioner(J, nblocks, gpu ? ExaPF.CUDADevice() : ExaPF.CPU())

    # Solve with basis vectors
    B = zeros(n, nrhs)
    for i in 1:nrhs
        B[i, i] = 1.0
    end

    if gpu
        # Transfer data to the GPU
        gJ = CuSparseMatrixCSR(J)
        gB = CuArray{Float64, 2}(B)
        @time LS.update(gJ, precond)
        (xsol, stat) = @time BPF.block_bicgstab(
            gJ, gB; M=precond.P, itmax=500, atol=1e-4
        )
    else
        @time LS.update(J, precond)
        (xsol, stat) = @time BPF.block_bicgstab(
            J, B; M=precond.P, itmax=500, atol=1e-4
        )
    end
    println(stat.status)
    println(length(stat.residuals))
    println(stat.residuals[end])
    return xsol, stat.residuals
end

datafile = joinpath(dirname(@__FILE__), "..", "data", "case300.txt")
J = mmread(datafile)

xsol, residuals = blk_bicgstab(J; nrhs=32, gpu=false)
