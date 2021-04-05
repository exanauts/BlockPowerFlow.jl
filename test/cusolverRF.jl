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
        CUDA.@time ldiv!(x, rflu, gb)

        # Test that first RHS is the same as the one computed
        # using UMFPACK
        res = Array(x)
        @test isapprox(res, solution)

        # Test update
        gJ.nzVal .*= 2.0
        CUSOLVERRF.update!(rflu, gJ)
        CUDA.@time ldiv!(x, rflu, gb)
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
        rflu = CUSOLVERRF.CusolverRfLUBatch(gJ, nbatch)

        # With vectors
        gB = Vector{CuVector{Float64}}(undef, nbatch)
        gX = Vector{CuVector{Float64}}(undef, nbatch)
        for i in 1:nbatch
            gB[i] = i .* gb
            gX[i] = copy(x)
        end
        CUDA.@time CUSOLVERRF.solve!(gX, rflu, gB)
        for i in 1:nbatch
            res = Array(gX[i])
            α = 1.0 * i
            @test isapprox(res, α * (J \ b[1:n]))
        end

        CUSOLVERRF.update!(rflu, gJ)

        # With matrices
        gBm = CuMatrix{Float64}(undef, n, nbatch)
        gXm = CuMatrix{Float64}(undef, n, nbatch)
        for i in 1:nbatch
            gBm[:, i] .= i .* gb
        end
        CUDA.@time ldiv!(gXm, rflu, gBm)
        for i in 1:nbatch
            res = Array(gXm[:, i])
            α = 1.0 * i
            @test isapprox(res, α * (J \ b[1:n]))
        end
        #
        CUSOLVERRF.cudestroy!(rflu)
    end
end

