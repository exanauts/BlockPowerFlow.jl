using Test
using LinearAlgebra
using MatrixMarket

using BlockPowerFlow

using CUDA
using CUDA.CUSPARSE

@testset "CusolverRF" begin
    datafile = joinpath(dirname(@__FILE__), "..", "data", "case300.txt")
    J = mmread(datafile)
    n, m = size(J)
    gJ = CuSparseMatrixCSR(J)
    b = rand(n)
    gb = CuVector{Float64}(b)
    x = CUDA.zeros(Float64, m)

    # Compute solution with UMFPACK
    solution = J \ b

    @testset "LU factorization" begin
        rflu = CUSOLVERRF.CusolverRfLU(gJ; ordering=:AMD)
        CUDA.@time CUSOLVERRF.solve!(x, rflu, gb)

        # Test that first RHS is the same as the one computed
        # using UMFPACK
        res = Array(x)
        @test isapprox(res, solution)

        # Test update
        gJ.nzVal .*= 2.0
        CUSOLVERRF.update!(rflu, gJ)
        CUDA.@time CUSOLVERRF.solve!(x, rflu, gb)
        res = Array(x)
        @test isapprox(res, 0.5 * solution)

        # Free
        CUSOLVERRF.cudestroy!(rflu)
    end

    # Reset matrix
    gJ = CuSparseMatrixCSR(J)
    @testset "LU factorization (batch)" begin
        # Test batch mode
        nbatch = 32
        gB = Vector{CuVector{Float64}}(undef, nbatch)
        gX = Vector{CuVector{Float64}}(undef, nbatch)
        for i in 1:nbatch
            gB[i] = gb
            gX[i] = x
        end
        rflu = CUSOLVERRF.CusolverRfLUBatch(gJ, nbatch)
        CUDA.@time CUSOLVERRF.solve!(gX, rflu, gB)
        for i in 1:nbatch
            res = Array(gX[i])
            α = 1.0
            @test isapprox(res, α * (J \ b[1:n]))
        end
        CUSOLVERRF.cudestroy!(rflu)
    end
end

