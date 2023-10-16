using MatrixMarket

using AMDGPU
using AMDGPU.rocSPARSE, SparseArrays
import KernelAbstractions
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
# name = "case9241pegase"
# name = "case13659pegase"
# name = "case_ACTIVSg10k"
# name = "case_ACTIVSg25k"
# name = "case_ACTIVSg70k"

names = ["case9", "case57", "case118", "case300", "case1354pegase",
         "case_ACTIVSg2000", "case2869pegase", "case9241pegase",
         "case13659pegase", "case_ACTIVSg10k", "case_ACTIVSg25k", "case_ACTIVSg70k"]

if AMDGPU.functional()
    for name in names
        path_A = joinpath(@__DIR__, "dump", "Gx_$(name).mtx")
        path_B = joinpath(@__DIR__, "dump", "Gu_$(name).mtx")

        # Storage on CPU
        A = mmread(path_A)
        B = mmread(path_B)
        B = Matrix(B)

        # Storage on GPU
        A = ROCSparseMatrixCSR(A)
        B = ROCMatrix(B)

        # Create the ILU(0) preconditioner
        P = BPF.ilu0(A)

        # Parameters
        m,n = size(A)
        s,p = size(B)
        verbose = 1
        rtol = 0.0
        atol = 1e-10
        memory = 1
        itmax = 50

        println("Problem $name of size ($m,$n) with $p right-hand sides.")

        # Solve the linear system with a right preconditioner ILU(0)
        println("Problem $name with a preconditioner")
        X_gmres, stats = BPF.block_gmres(A, B; N=P, ldiv=true, verbose, atol, rtol, memory, itmax)
        RNorm = norm(B - A * X_gmres)
        println("‖Rₖ‖: ", RNorm)

        # Solve the linear system without a preconditioner
        println("Problem $name without a preconditioner")
        X_gmres, stats = BPF.block_gmres(A, B; verbose, atol, rtol, memory, itmax)
        RNorm = norm(B - A * X_gmres)
        println("‖Rₖ‖: ", RNorm)
    end
end
