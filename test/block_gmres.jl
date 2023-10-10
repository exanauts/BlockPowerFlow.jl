using Test
using LinearAlgebra
using MatrixMarket

using BlockPowerFlow

# using CUDA
# using CUDA.CUSPARSE
import ExaPF
import ExaPF: LinearSolvers

const LS = LinearSolvers
const BPF = BlockPowerFlow

function blk_gmres(J; nrhs=10, nblocks=8, gpu=false)
    n = size(J, 1)
    # Init preconditioner
    device = gpu ? ExaPF.CUDADevice() : ExaPF.CPU()
    precond = LS.BlockJacobiPreconditioner(J, nblocks, device)

    # Solve with basis vectors
    B = zeros(n, nrhs)
    for i in 1:nrhs
        B[i, i] = 1.0
    end

    if gpu
        # Transfer data to the GPU
        gJ = CuSparseMatrixCSR(J)
        gB = CuArray{Float64, 2}(B)
        @time LS.update(precond, gJ, device)
        (xsol, stat) = @time BPF.block_gmres(
            gJ, gB; M=precond, history=true, itmax=500, atol=1e-4
        )
    else
        @time LS.update(precond, J, device)
        (xsol, stat) = @time BPF.block_gmres(
            J, B; M=precond, history=true, itmax=500, atol=1e-4
        )
    end
    println(stat.status)
    println(stat.niter)
    gpu ? println(norm(gB - gJ * xsol)) : println(norm(B - J * xsol))
    return xsol, stat.residuals
end

datafile = joinpath(dirname(@__FILE__), "..", "data", "case300.txt")
J = mmread(datafile)

xsol, residuals = blk_gmres(J; nrhs=32, gpu=false)
