
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
    x = CUDA.zeros(Float64, batch * m)

    rflu = CUSOLVERRF.CusolverRfLU(gJ)
    CUSOLVERRF.solve!(x, rflu, gb)

    # Test that first RHS is the same as the one computed
    # using UMFPACK
    res = Array(x[1:n])
    @test isapprox(res, J \ b[1:n])
end
