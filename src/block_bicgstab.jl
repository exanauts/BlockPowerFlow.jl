
# An implementation of block BICGSTAB for the solution of unsymmetric and square consistent linear system AX = B.
#
# This method is described in
#
# A. El Guennouni and K. Jbilou and H. Sadok
# A block version of BiCGSTAB for linear systems with multiple right-hand sides.
# Electronic Transactions on Numerical Analysis, 16, pp. 129--142, 2003.
#
# Alexis Montoison, <alexis.montoison@polymtl.ca>
# Montréal, December 2020.


mutable struct SimpleStats{T}
  solved :: Bool
  inconsistent :: Bool
  residuals :: Vector{T}
  Aresiduals :: Vector{T}
  status :: String
end

display2(iter, verbose) = (verbose > 0) && (mod(iter, verbose) == 0)

# export block_bicgstab

"""
    (X, stats) = block_bicgstab(A :: AbstractMatrix{T}, B :: AbstractMatrix{T}; C :: AbstractMatrix{T}=B,
                                M=opEye(), N=opEye(), atol :: T=√eps(T), rtol :: T=√eps(T),
                                itmax :: Int=0, verbose :: Int=0) where T <: AbstractFloat
"""
function block_bicgstab(A :: AbstractMatrix{T}, B :: AbstractMatrix{T}; C :: AbstractMatrix{T}=B,
                        M=opEye(), N=opEye(), atol :: T=√eps(T), rtol :: T=√eps(T),
                        itmax :: Int=0, verbose :: Int=0) where T <: AbstractFloat

  n, m = size(A)
  p, s = size(B)
  m == n || error("System must be square")
  m == p || error("Inconsistent problem size")
  (verbose > 0) && @printf("BLOCK-BICGSTAB: system of size %d with %d right-hand sides\n", n, s)

  # Check M == Iₘ and N == Iₙ
  MisI = isa(M, opEye)
  NisI = isa(N, opEye)

  # Check type consistency
  eltype(A) == T || error("eltype(A) ≠ $T")
  MisI || (eltype(M) == T) || error("eltype(M) ≠ $T")
  NisI || (eltype(N) == T) || error("eltype(N) ≠ $T")

  # Determine the storage type of B
  S = typeof(B)

  # Preconditionned matrix
  MA = M * A
  # Set up workspace.
  X = zeros(T, n, s)  # Xₖ
  R = similar(B)      # Rₖ = B - A * X
  mul!(R, M, B)
  mul!(R, MA, X, -1.0, 1.0)
  P = copy(R)         # Pₖ

  # Allocate working memory
  V = zeros(T, n, s)
  Q = zeros(T, n, s)
  K = zeros(T, s, s)
  α = zeros(T, s, s)
  β = zeros(T, s, s)

  # Compute residual norm ‖R₀‖₂.
  rNorm = norm(R)

  iter = 0
  itmax == 0 && (itmax = 2*n)

  rNorms = [rNorm;]
  ε = atol + rtol * rNorm
  (verbose > 0) && @printf("%5s  %7s\n", "k", "‖Rₖ‖")
  display2(iter, verbose) && @printf("%5d  %7.1e\n", iter, rNorm)

  # Stopping criterion.
  solved = rNorm ≤ ε
  tired = iter ≥ itmax
  breakdown = false
  status = "unknown"

  Cᵀ = copy(R)'

  while !(solved || tired)
    mul!(V, MA, P)              # V = M * A * P
    mul!(K, Cᵀ, V)              # K = C' * V

    # Dense factorization (inplace)
    K_factorized = lu!(K)

    mul!(α, Cᵀ, R)              # α = K \ (C' * R)
    ldiv!(K_factorized, α)

    mul!(R, V, α, -1.0, 1.0)    # R = R - V * α
    mul!(Q, MA, R)              # Q = M * A * R
    ω = norm(dot(Q, R)) / norm(dot(Q, Q))

    X .+= ω .* R                # X = X + ω * R
    mul!(X, P, α, 1.0, 1.0)     # X = X + P * α

    R .-= ω .* Q                # R = R - ω * Q

    mul!(β, Cᵀ, Q)              # β = K \ (C' * Q)
    ldiv!(K_factorized, β)

    # P = R + (P - ω * V) * β
    Q .= P .- ω .* V            # Q = P - ω * V
    P .= R                      # P = R
    mul!(P, Q, β, -1.0, 1.0)    # P = P - Q * β

    iter = iter + 1
    rNorm = norm(R)
    push!(rNorms, rNorm)

    # Update stopping criterion.
    solved = rNorm ≤ ε
    tired = iter ≥ itmax
    display2(iter, verbose) && @printf("%5d  %7.1e\n", iter, rNorm)
  end
  (verbose > 0) && @printf("\n")

  status = tired ? "maximum number of iterations exceeded" : "solution good enough given atol and rtol"
  stats = SimpleStats(solved, false, rNorms, T[], status)
  return (X, stats)
end

# We split the implementation on the GPU, as we are calling
# CUSOLVER directly in this case.
function block_bicgstab(A :: CuSparseMatrixCSR{T}, B :: CuMatrix{T}; C :: AbstractMatrix{T}=B,
                        M=opEye(), N=opEye(), atol :: T=√eps(T), rtol :: T=√eps(T),
                        itmax :: Int=0, verbose :: Int=0) where T <: AbstractFloat

  n, m = size(A)
  p, s = size(B)
  m == n || error("System must be square")
  m == p || error("Inconsistent problem size")
  (verbose > 0) && @printf("BLOCK-BICGSTAB: system of size %d with %d right-hand sides\n", n, s)

  # Check M == Iₘ and N == Iₙ
  MisI = isa(M, opEye)
  NisI = isa(N, opEye)

  # Check type consistency
  eltype(A) == T || error("eltype(A) ≠ $T")
  MisI || (eltype(M) == T) || error("eltype(M) ≠ $T")
  NisI || (eltype(N) == T) || error("eltype(N) ≠ $T")

  # Determine the storage type of B
  S = typeof(B)

  # Preconditionned matrix
  Ax = CUDA.zeros(T, n, s)
  # Set up workspace.
  X = CUDA.zeros(T, n, s)  # Xₖ
  R = similar(B)      # Rₖ = B - A * X
  mul!(R, M, B)
  P = copy(R)         # Pₖ

  # Allocate working memory
  V = CUDA.zeros(T, n, s)
  Q = CUDA.zeros(T, n, s)
  K = CUDA.zeros(T, s, s)
  α = CUDA.zeros(T, s, s)
  β = CUDA.zeros(T, s, s)

  # Compute residual norm ‖R₀‖₂.
  rNorm = norm(R)

  iter = 0
  itmax == 0 && (itmax = 2*n)

  rNorms = [rNorm;]
  ε = atol + rtol * rNorm
  (verbose > 0) && @printf("%5s  %7s\n", "k", "‖Rₖ‖")
  display2(iter, verbose) && @printf("%5d  %7.1e\n", iter, rNorm)

  # Stopping criterion.
  solved = rNorm ≤ ε
  tired = iter ≥ itmax
  breakdown = false
  status = "unknown"

  Cᵀ = copy(R)'

  while !(solved || tired)
    mul!(Ax, A, P)
    mul!(V, M, Ax)              # V = M * A * P
    mul!(K, Cᵀ, V)              # K = C' * V

    # Dense factorization (inplace)
    K_factorized, perm = CUSOLVER.getrf!(K)

    mul!(α, Cᵀ, R)              # α = K \ (C' * R)
    α = CUSOLVER.getrs!('N', K_factorized, perm, α)

    mul!(R, V, α, -1.0, 1.0)    # R = R - V * α
    mul!(Ax, A, R)
    mul!(Q, M, Ax)              # Q = M * A * R
    ω = norm(dot(Q, R)) / norm(dot(Q, Q))

    X .+= ω .* R                # X = X + ω * R
    mul!(X, P, α, 1.0, 1.0)     # X = X + P * α

    R .-= ω .* Q                # R = R - ω * Q

    mul!(β, Cᵀ, Q)              # β = K \ (C' * Q)
    β = CUSOLVER.getrs!('N', K_factorized, perm, β)

    # P = R + (P - ω * V) * β
    Q .= P .- ω .* V            # Q = P - ω * V
    P .= R                      # P = R
    mul!(P, Q, β, -1.0, 1.0)    # P = P - Q * β

    iter = iter + 1
    rNorm = norm(R)
    push!(rNorms, rNorm)

    # Update stopping criterion.
    solved = rNorm ≤ ε
    tired = iter ≥ itmax
    display2(iter, verbose) && @printf("%5d  %7.1e\n", iter, rNorm)
  end
  (verbose > 0) && @printf("\n")

  status = tired ? "maximum number of iterations exceeded" : "solution good enough given atol and rtol"
  stats = SimpleStats(solved, false, rNorms, T[], status)
  return (X, stats)
end

