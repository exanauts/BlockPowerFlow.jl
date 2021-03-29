
@testset "LU factorization (naive API)" begin
    datafile = joinpath(dirname(@__FILE__), "..", "data", "case300.txt")
    J = mmread(datafile)
    n, m = size(J)
    batch = 1
    gJ = CuSparseMatrixCSR(J)
    b = rand(batch * n)
    gb = CuVector{Float64}(b)
    x = CUDA.zeros(Float64, batch * m)

    CUDA.@time CUSOLVERRF.batch_lu!(gJ, gb, x, 1e-8, Cint(batch), 'O')

    # Test that first RHS is the same as the one computed
    # using UMFPACK
    res = Array(x[1:n])
    @test isapprox(res, J \ b[1:n])

end

@testset "LU factorization (new API)" begin
    datafile = joinpath(dirname(@__FILE__), "..", "data", "case300.txt")
    J = mmread(datafile)
    n, m = size(J)
    batch = 1
    gJ = CuSparseMatrixCSR(J)
    b = rand(batch * n)
    gb = CuVector{Float64}(b)
    x = CUDA.zeros(Float64, m)

    rflu = CUSOLVERRF.CusolverRfLU(gJ; ordering=:AMD)
    CUDA.@time CUSOLVERRF.solve!(x, rflu, gb)

    # Test that first RHS is the same as the one computed
    # using UMFPACK
    res = Array(x[1:n])
    @test isapprox(res, J \ b[1:n])

    # Test update
    gJ.nzVal .*= 2.0
    CUSOLVERRF.update!(rflu, gJ)
    CUDA.@time CUSOLVERRF.solve!(x, rflu, gb)
    res = Array(x[1:n])
    @test isapprox(res, 0.5 * (J \ b[1:n]))

    # Free
    CUSOLVERRF.cudestroy!(rflu)

    gJ.nzVal .*= 0.5
    # Test batch mode
    nbatch = 32
    gB = Vector{CuVector{Float64}}(undef, nbatch)
    gX = Vector{CuVector{Float64}}(undef, nbatch)
    for i in 1:nbatch
        gB[i] = gb
        gX[i] = x
    end
    rflu = CUSOLVERRF.CusolverRfLU(gJ; batchsize=nbatch)
    CUDA.@time CUSOLVERRF.solve_batch!(gX, rflu, gB)
    for i in 1:nbatch
        res = Array(gX[i])
        @test isapprox(res, (J \ b[1:n]))
    end
    CUSOLVERRF.cudestroy!(rflu)
end

