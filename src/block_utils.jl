"""
    FloatOrComplex{T}
Union type of `T` and `Complex{T}` where T is an `AbstractFloat`.
"""
const FloatOrComplex{T} = Union{T, Complex{T}} where T <: AbstractFloat

# Reduced QR factorization with Householder reflections:
# Q, R = householder(A)
#
# Input :
# A an n-by-k matrix, n ≥ k
#
# Output :
# Q an n-by-k orthonormal matrix: QᴴQ = Iₖ
# R an k-by-k upper triangular matrix: QR = A
function householder(A::AbstractMatrix{FC}; compact::Bool=false) where FC <: FloatOrComplex
    n, k = size(A)
    Q = copy(A)
    τ = zeros(FC, k)
    R = zeros(FC, k, k)
    householder!(Q, R, τ; compact)
end

function householder!(Q::AbstractMatrix{FC}, R::AbstractMatrix{FC}, τ::AbstractVector{FC}; compact::Bool=false) where FC <: FloatOrComplex
    n, k = size(Q)
    R .= zero(FC)
    LAPACK.geqrf!(Q, τ)
    for i = 1:k
        for j = i:k
            R[i,j] = Q[i,j]
        end
    end
    !compact && LAPACK.orgqr!(Q, τ)
    return Q, R
end

kdisplay(iter::Int, verbose::Int) = (verbose > 0) && (mod(iter, verbose) == 0)

ktimer(start_time::UInt64) = (time_ns() - start_time) / 1e9

"Abstract type for statistics returned by a solver"
abstract type KrylovStats{T} end

mutable struct BlockGMRESStats{T} <: KrylovStats{T}
    niter     :: Int
    solved    :: Bool
    residuals :: Vector{T}
    timer     :: Float64
    status    :: String
end

function reset!(stats :: BlockGMRESStats)
    empty!(stats.residuals)
end

"""
    M = vector_to_matrix(S)

Return the dense matrix storage type `M` related to the dense vector storage type `S`.
"""
function vector_to_matrix(::Type{S}) where S <: DenseVector
    T = hasproperty(S, :body) ? S.body : S
    par = T.parameters
    npar = length(par)
    (2 ≤ npar ≤ 3) || error("Type $S is not supported.")
    if npar == 2
      M = T.name.wrapper{par[1], 2}
    else
      M = T.name.wrapper{par[1], 2, par[3]}
    end
    return M
end

"""
    S = matrix_to_vector(M)

Return the dense vector storage type `S` related to the dense matrix storage type `M`.
"""
function matrix_to_vector(::Type{M}) where M <: DenseMatrix
    T = hasproperty(M, :body) ? M.body : M
    par = T.parameters
    npar = length(par)
    (2 ≤ npar ≤ 3) || error("Type $M is not supported.")
    if npar == 2
        S = T.name.wrapper{par[1], 1}
    else
        S = T.name.wrapper{par[1], 1, par[3]}
    end
    return S
end

allocate_if(bool, solver, v, S, n) = bool && isempty(solver.:($v)::S) && (solver.:($v)::S = S(undef, n))
allocate_if(bool, solver, v, S, m, n) = bool && isempty(solver.:($v)::S) && (solver.:($v)::S = S(undef, m, n))

mulorldiv!(y, P, x, ldiv::Bool) = ldiv ? ldiv!(y, P, x) : mul!(y, P, x)
