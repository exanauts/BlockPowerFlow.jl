
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

export BlockBicgstabSolver, block_bicgstab, block_bicgstab!

mutable struct SimpleStats{T}
  niter     :: Int
  solved    :: Bool
  residuals :: Vector{T}
  timer     :: Float64
  status    :: String
end

abstract type KrylovSolver{T,S} end

struct BlockBicgstabSolver{T,S} <: KrylovSolver{T,S}
  X     :: S
  Y     :: S
  P     :: S
  R     :: S
  V     :: S
  Q     :: S
  U     :: S
  K     :: S
  D     :: S
  α     :: S
  β     :: S
  stats :: SimpleStats{T}
end

function BlockBicgstabSolver(A, B)
  n, m = size(A)
  s, p = size(B)
  S = typeof(B)
  T = eltype(B)
  X = S(undef, n, p)
  Y = S(undef, n, p)
  P = S(undef, n, p)
  R = S(undef, n, p)
  V = S(undef, n, p)
  Q = S(undef, n, p)
  U = S(undef, n, p)
  K = S(undef, p, p)
  D = S(undef, p, p)
  α = S(undef, p, p)
  β = S(undef, p, p)
  stats = SimpleStats(0, false, T[], 0.0, "unknown")
  return BlockBicgstabSolver{T,S}(X, Y, P, R, V, Q, U, K, D, α, β, stats)
end

"""
    (X, stats) = block_bicgstab(A :: AbstractMatrix{T}, B :: AbstractMatrix{T}; C :: AbstractMatrix{T}=B,
                                N=opEye(), atol :: T=√eps(T), rtol :: T=√eps(T),
                                itmax :: Int=0, verbose :: Int=0, history :: Bool=false) where T <: AbstractFloat
"""
function block_bicgstab(A :: AbstractMatrix{T}, B :: AbstractMatrix{T}; kwargs...) where T <: AbstractFloat
  solver = BlockBicgstabSolver(A, B)
  block_bicgstab!(solver, A, B; kwargs...)
end

function block_bicgstab!(solver :: BlockBicgstabSolver{T,Matrix{T}},
                         A :: AbstractMatrix{T}, B :: AbstractMatrix{T}; C :: AbstractMatrix{T}=B,
                         N=I, atol :: T=√eps(T), rtol :: T=√eps(T),
                         itmax :: Int=0, timemax::Float64=Inf, verbose :: Int=0, history :: Bool=false) where T <: AbstractFloat

  # Timer
  start_time = time_ns()
  timemax_ns = 1e9 * timemax

  n, m = size(A)
  s, p = size(B)
  m == n || error("System must be square")
  n == s || error("Inconsistent problem size")
  (verbose > 0) && @printf("BLOCK-BICGSTAB: system of size %d with %d right-hand sides\n", n, p)

  # Check N == Iₙ
  NisI = (N === I)

  # Check type consistency
  eltype(A) == T || error("eltype(A) ≠ $T")
  NisI || (eltype(N) == T) || error("eltype(N) ≠ $T")

  # Set up workspace.
  X, P, R, V, Q, U, K, D, α, β, stats = solver.X, solver.P, solver.R, solver.V, solver.Q, solver.U, solver.K, solver.D, solver.α, solver.β, solver.stats
  Cᵀ = C'

  X .= zero(T)
  P .= B
  R .= B
  V .= zero(T)
  Q .= zero(T)
  U .= zero(T)
  K .= zero(T)
  D .= zero(T)
  α .= zero(T)
  β .= zero(T)

  # Right preconditioner
  if !NisI
    Y = solver.Y
    Y .= zero(T)
  end

  # Compute residual norm ‖R₀‖₂.
  rNorm = norm(R)

  iter = 0
  itmax == 0 && (itmax = 2*n)

  rNorms = stats.residuals
  !history && !isempty(rNorms) && (rNorms = T[])
  history && push!(rNorms, rNorm)
  ε = atol + rtol * rNorm
  (verbose > 0) && @printf("%5s  %7s  %5s\n", "k", "‖Rₖ‖", "timer")
  kdisplay(iter, verbose) && @printf("%5d  %7.1e  %.2fs\n", iter, rNorm, ktimer(start_time))

  # Stopping criterion.
  solved = rNorm ≤ ε
  tired = iter ≥ itmax
  overtimed = false
  status = "unknown"

  while !(solved || tired || overtimed)

    NisI ? Y = P : mul!(Y, N, P)           # Yₖ = N⁻¹Pₖ
    mul!(V, A, Y)                          # Vₖ = AYₖ
    mul!(K, Cᵀ, V)                         # Kₖ = CᵀVₖ
    mul!(D, Cᵀ, R)                         # Dₖ = CᵀRₖ
    K_factorized = lu!(K)
    α .= D; ldiv!(K_factorized, α)         # Kₖαₖ = Dₖ, αₖ is a p×p matrix
    U .= R
    mul!(U, V, α, -one(T), one(T))         # Uₖ = Rₖ - Vₖαₖ
    mul!(R, Y, α)                          # Rₐᵤₓ = N⁻¹Pₖαₖ, temporary storage
    @. X = X + R                           # Xₐᵤₓ = Xₖ + Rₐᵤₓ
    NisI ? Y = U : mul!(Y, N, U)           # Yₖ = N⁻¹Uₖ
    mul!(Q, A, Y)                          # Qₖ = AYₖ
    ω = norm(dot(Q, U)) / norm(dot(Q, Q))  # ωₖ = ‖⟨Qₖ,Uₖ⟩‖ / ‖⟨Qₖ,Qₖ⟩‖, ‖•‖ is the Frobenius norm
    @. X = X + ω * Y                       # Xₖ₊₁ = Xₐᵤₓ + ωₖN⁻¹Uₖ
    @. R = U - ω * Q                       # Rₖ₊₁ = Uₖ - ωₖQₖ
    D = mul!(D, Cᵀ, Q)                     # Dₖ = CᵀQₖ
    β .= D; ldiv!(K_factorized, β)         # Kₖβₖ = Dₖ, βₖ is a p×p matrix
    @. U = P - ω * V                       # Uₐᵤₓ = Pₖ - ωₖVₖ, temporary storage
    mul!(P, U, β)                          # Pₐᵤₓ = Uₖβₖ
    @. P = R - P                           # Pₖ₊₁ = Rₖ₊₁ - Pₐᵤₓ

    iter = iter + 1
    rNorm = norm(R)
    history && push!(rNorms, rNorm)

    # Update stopping criterion.
    solved = rNorm ≤ ε
    tired = iter ≥ itmax
    timer = time_ns() - start_time
    overtimed = timer > timemax_ns
    kdisplay(iter, verbose) && @printf("%5d  %7.1e  %.2fs\n", iter, rNorm, ktimer(start_time))
  end
  (verbose > 0) && @printf("\n")

  # Termination status
  tired     && (status = "maximum number of iterations exceeded")
  solved    && (status = "solution good enough given atol and rtol")
  overtimed && (status = "time limit exceeded")

  stats.niter = iter
  stats.solved = solved
  stats.timer = ktimer(start_time)
  stats.status = status
  return (X, stats)
end

# We split the implementation on the GPU, as we are calling
# CUSOLVER directly in this case.
function block_bicgstab(A :: CuSparseMatrixCSR{T}, B :: CuMatrix{T}; kwargs...) where T <: AbstractFloat
  solver = BlockBicgstabSolver(A, B)
  block_bicgstab!(solver, A, B; kwargs...)
end

function block_bicgstab!(solver :: BlockBicgstabSolver{T,<:CuMatrix{T}},
                         A :: CuSparseMatrixCSR{T}, B :: CuMatrix{T}; C :: AbstractMatrix{T}=B,
                         N=I, atol :: T=√eps(T), rtol :: T=√eps(T),
                         itmax :: Int=0, timemax::Float64=Inf, verbose :: Int=0, history :: Bool=false) where T <: AbstractFloat

  # Timer
  start_time = time_ns()
  timemax_ns = 1e9 * timemax

  n, m = size(A)
  s, p = size(B)
  m == n || error("System must be square")
  n == s || error("Inconsistent problem size")
  (verbose > 0) && @printf("BLOCK-BICGSTAB: system of size %d with %d right-hand sides\n", n, p)

  # Check N == Iₙ
  NisI = (N === I)

  # Check type consistency
  eltype(A) == T || error("eltype(A) ≠ $T")
  NisI || (eltype(N) == T) || error("eltype(N) ≠ $T")

  # Set up workspace.
  X, P, R, V, Q, U, K, D, α, β, stats = solver.X, solver.P, solver.R, solver.V, solver.Q, solver.U, solver.K, solver.D, solver.α, solver.β, solver.stats
  Cᵀ = C'

  X .= zero(T)
  P .= B
  R .= B
  V .= zero(T)
  Q .= zero(T)
  U .= zero(T)
  K .= zero(T)
  D .= zero(T)
  α .= zero(T)
  β .= zero(T)

  # Right preconditioner
  if !NisI
    Y = solver.Y
    Y .= zero(T)
  end

  # Compute residual norm ‖R₀‖₂.
  rNorm = norm(R)

  iter = 0
  itmax == 0 && (itmax = 2*n)

  rNorms = stats.residuals
  !history && !isempty(rNorms) && (rNorms = T[])
  history && push!(rNorms, rNorm)
  ε = atol + rtol * rNorm
  (verbose > 0) && @printf("%5s  %7s  %5s\n", "k", "‖Rₖ‖", "timer")
  kdisplay(iter, verbose) && @printf("%5d  %7.1e  %.2fs\n", iter, rNorm, ktimer(start_time))

  # Stopping criterion.
  solved = rNorm ≤ ε
  tired = iter ≥ itmax
  overtimed = false
  status = "unknown"

  while !(solved || tired || overtimed)

    NisI ? Y = P : mul!(Y, N, P)           # Yₖ = N⁻¹Pₖ
    mul!(V, A, Y)                          # Vₖ = AYₖ
    mul!(K, Cᵀ, V)                         # Kₖ = CᵀVₖ
    mul!(D, Cᵀ, R)                         # Dₖ = CᵀRₖ
    α .= D                                 # Kₖαₖ = Dₖ, αₖ is a p×p matrix
    K_factorized, τ = CUSOLVER.geqrf!(K)
    CUSOLVER.ormqr!('L', 'T', K_factorized, τ, α)
    CUBLAS.trsm!('L', 'U', 'N', 'N', one(T), K_factorized, α)
    U .= R
    mul!(U, V, α, -one(T), one(T))         # Uₖ = Rₖ - Vₖαₖ
    mul!(R, Y, α)                          # Rₐᵤₓ = N⁻¹Pₖαₖ, temporary storage
    X .+= R                                # Xₐᵤₓ = Xₖ + Rₐᵤₓ
    NisI ? Y = U : mul!(Y, N, U)           # Yₖ = N⁻¹Uₖ
    mul!(Q, A, Y)                          # Qₖ = AYₖ
    ω = norm(dot(Q, U)) / norm(dot(Q, Q))  # ωₖ = ‖⟨Qₖ,Uₖ⟩‖ / ‖⟨Qₖ,Qₖ⟩‖, ‖•‖ is the Frobenius norm
    X .+= ω .* Y                           # Xₖ₊₁ = Xₐᵤₓ + ωₖN⁻¹Uₖ
    R .= U .- ω .* Q                       # Rₖ₊₁ = Uₖ - ωₖQₖ
    D = mul!(D, Cᵀ, Q)                     # Dₖ = CᵀQₖ
    β .= D                                 # Kₖβₖ = Dₖ, βₖ is a p×p matrix
    CUSOLVER.ormqr!('L', 'T', K_factorized, τ, β)
    CUBLAS.trsm!('L', 'U', 'N', 'N', one(T), K_factorized, β)
    U .= P .- ω .* V                       # Uₐᵤₓ = Pₖ - ωₖVₖ, temporary storage
    mul!(P, U, β)                          # Pₐᵤₓ = Uₖβₖ
    P .= R .- P                            # Pₖ₊₁ = Rₖ₊₁ - Pₐᵤₓ

    iter = iter + 1
    rNorm = norm(R)
    history && push!(rNorms, rNorm)

    # Update stopping criterion.
    solved = rNorm ≤ ε
    tired = iter ≥ itmax
    timer = time_ns() - start_time
    overtimed = timer > timemax_ns
    kdisplay(iter, verbose) && @printf("%5d  %7.1e  %.2fs\n", iter, rNorm, ktimer(start_time))
  end
  (verbose > 0) && @printf("\n")

  # Termination status
  tired     && (status = "maximum number of iterations exceeded")
  solved    && (status = "solution good enough given atol and rtol")
  overtimed && (status = "time limit exceeded")

  stats.niter = iter
  stats.solved = solved
  stats.timer = ktimer(start_time)
  stats.status = status
  return (X, stats)
end
