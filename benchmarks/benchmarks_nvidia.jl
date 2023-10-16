using MatrixMarket

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

path_A = joinpath(@__DIR__, "dump", "Gx_$(name).mtx")
path_B = joinpath(@__DIR__, "dump", "Gu_$(name).mtx")

A = mmread(path_A) 
B = mmread(path_B)
B = Matrix(B)

verbose = 1
m, n = size(A)
itmax = 5
gpu = false

nblocks = div(n, 32)
device = gpu ? ExaPF.CUDABackend() : ExaPF.CPU()
P = LS.BlockJacobiPreconditioner(A, nblocks, device)
A = gpu ? CuSparseMatrixCSR(A) : A
B = gpu ? CuMatrix(B) : B
LS.update(P, A, device)

X_gmres, stats = BPF.block_gmres(A, B; N=P, memory=1, verbose, itmax, rtol=0.0, atol=1e-10)
RNorm = norm(B - A * X_gmres)
println("‖Rₖ‖ : ", RNorm)
