
# Reduced QR factorization with Householder reflections:
# Q, R = householder(A)
#
# Input :
# A an n-by-k matrix, n ≥ k
#
# Output :
# Q an n-by-k orthonormal matrix: QᴴQ = Iₖ
# R an k-by-k upper triangular matrix: QR = A
function householder(A::AbstractMatrix{FC}) where FC <: FloatOrComplex
  n, k = size(A)
  Q = copy(A)
  τ = zeros(FC, k)
  R = zeros(FC, k, k)
  householder!(Q, R, τ)
end

function householder!(Q::AbstractMatrix{FC}, R::AbstractMatrix{FC}, τ::AbstractVector{FC}) where FC <: FloatOrComplex
  n, k = size(Q)
  R .= zero(FC)
  @kgeqrf!(Q, τ)
  for i = 1:k
    for j = i:k
      R[i,j] = Q[i,j]
    end
  end
  @korgqr!(Q, τ)
  return Q, R
end

"""
    FloatOrComplex{T}
Union type of `T` and `Complex{T}` where T is an `AbstractFloat`.
"""
const FloatOrComplex{T} = Union{T, Complex{T}} where T <: AbstractFloat

kdisplay(iter::Int, verbose::Bool) = (verbose > 0) && (mod(iter, verbose) == 0)

ktimer(start_time::UInt64) = (time_ns() - start_time) / 1e9

"Abstract type for statistics returned by a solver"
abstract type KrylovStats{T} end

mutable struct BlockGmresStats{T} <: KrylovStats{T}
  niter     :: Int
  solved    :: Bool
  residuals :: Vector{T}
  timer     :: Float64
  status    :: String
end

function reset!(stats :: BlockGmresStats)
  empty!(stats.residuals)
end

abstract type KrylovSolver{T,FC,S} end
