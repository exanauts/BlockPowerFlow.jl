# An implementation of block-GMRES for the solution of the square linear system AX = B.
#
# Alexis Montoison, <alexis.montoison@polymtl.ca>
# Chicago, October 2023.

include("block_utils.jl")

export BlockGmresSolver, block_gmres, block_gmres!

"""
Type for storing the vectors required by the in-place version of BLOCK-GMRES.

The outer constructors

    solver = BlockGmresSolver(m, n, p, memory, S)
    solver = BlockGmresSolver(A, B, memory=5)

may be used in order to create these vectors.
`memory` is set to `div(n,p)` if the value given is larger than `div(n,p)`.
"""
mutable struct BlockGmresSolver{T,FC,S} <: KrylovSolver{T,FC,S}
  m          :: Int
  n          :: Int
  p          :: Int
  ΔX         :: S
  X          :: S
  W          :: S
  P          :: S
  Q          :: S
  V          :: Vector{S}
  Z          :: Vector{S}
  R          :: Vector{S}
  warm_start :: Bool
  stats      :: BlockGmresStats{T}
end

function BlockGmresSolver(m, n, p, memory, S)
  memory = min(div(n,p), memory)
  FC = eltype(S)
  T  = real(FC)
  ΔX = S(undef, 0, 0)
  X  = S(undef, n, p)
  W  = S(undef, n, p)
  P  = S(undef, 0)
  Q  = S(undef, 0)
  V  = S[S(undef, n, p) for i = 1 : memory]
  Z  = S[S(undef, p, p) for i = 1 : memory]
  R  = S[S(undef, p, p) for i = 1 : div(memory * (memory+1), 2)]
  stats = BlockGmresStats(0, false, T[], 0.0, "unknown")
  solver = BlockGmresSolver{T,FC,S}(m, n, p, ΔX, X, W, P, Q, V, Z, R, false, stats)
  return solver
end

function GmresSolver(A, B, memory=5)
  m, n = size(A)
  s, p = size(B)
  S = typeof(B)
  BlockGmresSolver(m, n, p, memory, S)
end

"""
    (X, stats) = block_gmres(A, B; X0::AbstractMatrix{FC}; memory::Int=20, M=I, N=I,
                             ldiv::Bool=false, restart::Bool=false, reorthogonalization::Bool=false,
                             atol::T = √eps(T), rtol::T=√eps(T), itmax::Int=0,
                             timemax::Float64=Inf, verbose::Int=0, history::Bool=false) where {T <: AbstractFloat, FC <: FloatOrComplex{T}}
"""
function block_gmres end

function block_gmres(A, B::AbstractMatrix{FC}, X0::AbstractMatrix{FC}; memory::Int=20, M=I, N=I,
                     ldiv::Bool=false, restart::Bool=false, reorthogonalization::Bool=false,
                     atol::T = √eps(T), rtol::T=√eps(T), itmax::Int=0,
                     timemax::Float64=Inf, verbose::Int=0, history::Bool=false) where {T <: AbstractFloat, FC <: FloatOrComplex{T}}

  start_time = time_ns()
  solver = blockGmresSolver(A, B, memory)
  warm_start!(solver, X0)
  elapsed_time = ktimer(start_time)
  timemax -= elapsed_time
  block_gmres!(solver, A, B; M, N, ldiv, restart, reorthogonalization, atol, rtol, itmax, timemax, verbose, history)
  solver.stats.timer += elapsed_time
  return solver.X, solver.stats
end

function block_gmres(A, B::AbstractMatrix{FC}; memory::Int=20, M=I, N=I,
                     ldiv::Bool=false, restart::Bool=false, reorthogonalization::Bool=false,
                     atol::T = √eps(T), rtol::T=√eps(T), itmax::Int=0,
                     timemax::Float64=Inf, verbose::Int=0, history::Bool=false) where {T <: AbstractFloat, FC <: FloatOrComplex{T}}

  start_time = time_ns()
  solver = BlockGmresSolver(A, B, memory)
  elapsed_time = ktimer(start_time)
  timemax -= elapsed_time
  block_gmres!(solver, A, B; M, N, ldiv, restart, reorthogonalization, atol, rtol, itmax, timemax, verbose, history)
  solver.stats.timer += elapsed_time
  return solver.X, solver.stats
end

function block_gmres!(solver :: GmresSolver{T,FC,S}, A::AbstractMatrix{FC}, X0::AbstractMatrix{FC}; M=I, N=I,
                      ldiv::Bool=false, restart::Bool=false, reorthogonalization::Bool=false,
                      atol::T = √eps(T), rtol::T=√eps(T), itmax::Int=0,
                      timemax::Float64=Inf, verbose::Int=0, history::Bool=false) where {T <: AbstractFloat, FC <: FloatOrComplex{T}, S <: AbstractMatrix{FC}}

  start_time = time_ns()
  warm_start!(solver, X0)
  elapsed_time = ktimer(start_time)
  timemax -= elapsed_time
  block_gmres!(solver, A, B; M, N, ldiv, restart, reorthogonalization, atol, rtol, itmax, timemax, verbose, history)
  solver.stats.timer += elapsed_time
  return solver
end

function block_gmres!(solver :: GmresSolver{T,FC,S}, A, B::AbstractMatrix{FC}; M=I, N=I,
                      ldiv::Bool=false, restart::Bool=false, reorthogonalization::Bool=false,
                      atol::T = √eps(T), rtol::T=√eps(T), itmax::Int=0,
                      timemax::Float64=Inf, verbose::Int=0, history::Bool=false) where {T <: AbstractFloat, FC <: FloatOrComplex{T}, S <: AbstractMatrix{FC}}

  # Timer
  start_time = time_ns()
  timemax_ns = 1e9 * timemax

  n, m = size(A)
  s, p = size(B)
  m == n || error("System must be square")
  n == s || error("Inconsistent problem size")
  (verbose > 0) && @printf("BLOCK-GMRES: system of size %d with %d right-hand sides\n", n, p)

  # Check M = Iₙ and N = Iₙ
  MisI = (M === I)
  NisI = (N === I)

  # Check type consistency
  eltype(A) == FC || @warn "eltype(A) ≠ $FC. This could lead to errors or additional allocations in operator-matrix products."
  ktypeof(B) <: S || error("ktypeof(B) is not a subtype of $S")

  # Set up workspace.
  allocate_if(!MisI  , solver, :Q , S, n, p)
  allocate_if(!NisI  , solver, :P , S, n, p)
  allocate_if(restart, solver, :ΔX, S, n, p)
  ΔX, X, W, V, Z = solver.ΔX, solver.X, solver.W, solver.V, solver.Z
  R, stats = solver.R, solver.stats
  warm_start = solver.warm_start
  RNorms = stats.residuals
  reset!(stats)
  Q  = MisI ? W : solver.Q
  R₀ = MisI ? W : solver.Q
  Xr = restart ? ΔX : X

  # FIX ME
  Ψtmp = Ztmp = S(undef, p, p)
  τ = S(undef, p)

  # Coefficients for mul!
  α = -one(FC)
  β = one(FC)
  γ = one(FC)

  # Initial solution X₀.
  X .= zero(FC)

  # Initial residual R₀.
  if warm_start
    mul!(W, A, Δx)
    W .= B .- W
    restart && (X .+= ΔX)
  else
    W .= B
  end
  MisI || mulorldiv!(R₀, M, W, ldiv)  # R₀ = M(B - AX₀)
  RNorm = norm(R₀)                    # ‖R₀‖_F   

  history && push!(RNorms, RNorm)
  ε = atol + rtol * RNorm

  mem = length(V)  # Memory
  npass = 0        # Number of pass

  iter = 0        # Cumulative number of iterations
  inner_iter = 0  # Number of iterations in a pass

  itmax == 0 && (itmax = 2*div(n,p))
  inner_itmax = itmax

  (verbose > 0) && @printf("%5s  %5s  %7s  %5s\n", "pass", "k", "‖Rₖ‖", "timer")
  kdisplay(iter, verbose) && @printf("%5d  %5d  %7.1e  %.2fs\n", npass, iter, rNorm, ktimer(start_time))

  # Stopping criterion
  solved = rNorm ≤ ε
  tired = iter ≥ itmax
  inner_tired = inner_iter ≥ inner_itmax
  status = "unknown"
  overtimed = false

  while !(solved || tired ||  overtimed)

    # Initialize workspace.
    nr = 0  # Number of blocks Ψᵢⱼ stored in Rₖ.
    for i = 1 : mem
      V[i] .= zero(FC)  # Orthogonal basis of Kₖ(MAN, MR₀).
    end
    for Ψ in R
      block .= zero(FC)  # Upper triangular matrix Rₖ.
    end
    for block in Z
      block .= zero(FC)  # Right-hand of the least squares problem min ‖Hₖ₊₁.ₖYₖ - ΓE₁‖₂.
    end

    if restart
      Xr .= zero(FC)  # Xr === ΔX when restart is set to true
      if npass ≥ 1
        mul!(W, A, X)
        W .= B .- W
        MisI || mulorldiv!(R₀, M, W, ldiv)
      end
    end
    
    # Initial Γ and V₁
    V[1] .= R₀
    householder!(V[1], Z[1], τ)

    npass = npass + 1
    solver.inner_iter = 0
    inner_tired = false

    while !(solved || inner_tired || overtimed)

      # Update iteration index
      solver.inner_iter = solver.inner_iter + 1
      inner_iter = solver.inner_iter

      # Update workspace if more storage is required and restart is set to false
      if !restart && (inner_iter > mem)
        for i = 1 : inner_iter
          push!(R, S(undef, p, p))
        end
      end

      # Continue the block-Arnoldi process.
      P = NisI ? V[inner_iter] : solver.P
      NisI || mulorldiv!(P, N, V[inner_iter], ldiv)  # P ← NVₖ
      mul!(W, A, P)                                  # W ← ANVₖ
      MisI || mulorldiv!(Q, M, W, ldiv)              # Q ← MANVₖ
      for i = 1 : inner_iter
        mul!(R[nr+i], V[i]', Q)       # Ψᵢₖ = Vᵢᴴ * Q
        mul!(Q, V[i], R[nr+i], α, β)  # Q = Q - Vᵢ * Ψᵢₖ
      end

      # Reorthogonalization of the block-Krylov basis.
      if reorthogonalization
        for i = 1 : inner_iter
          mul!(Ψtmp, V[i]', Q)       # Ψtmp = Vᵢᴴ * Q
          mul!(Q, V[i], Ψtmp, α, β)  # Q = Q - Vᵢ * Ψtmp
          R[nr+i] .+= Ψtmp
        end
      end

      # Vₖ₊₁ and Ψₖ₊₁.ₖ are stored in Q and Ztmp. 
      householder!(Q, Ztmp, τ)

      # Update the QR factorization of Hₖ₊₁.ₖ.
      # Apply previous Givens reflections Ωᵢ.
      # [cᵢ  sᵢ] [ r̄ᵢ.ₖ ] = [ rᵢ.ₖ ]
      # [s̄ᵢ -cᵢ] [rᵢ₊₁.ₖ]   [r̄ᵢ₊₁.ₖ]
      for i = 1 : inner_iter-1
        Rtmp      =      c[i]  * R[nr+i] + s[i] * R[nr+i+1]
        R[nr+i+1] = conj(s[i]) * R[nr+i] - c[i] * R[nr+i+1]
        R[nr+i]   = Rtmp
      end

      # Compute and apply current Givens reflection Ωₖ.
      # [cₖ  sₖ] [ r̄ₖ.ₖ ] = [rₖ.ₖ]
      # [s̄ₖ -cₖ] [hₖ₊₁.ₖ]   [ 0  ]
      (c[inner_iter], s[inner_iter], R[nr+inner_iter]) = sym_givens(R[nr+inner_iter], Hbis)

      # Update Zₖ = (Qₖ)ᴴΓE₁
      ζₖ₊₁          = conj(s[inner_iter]) * z[inner_iter]
      z[inner_iter] =      c[inner_iter]  * z[inner_iter]

      # Update residual norm estimate.
      # ‖ M(B - AXₖ) ‖_F = ‖Ztmp‖_F
      RNorm = norm(Ztmp)
      history && push!(RNorms, RNorm)

      # Update the number of coefficients in Rₖ
      nr = nr + inner_iter

      # Update stopping criterion.
      solved = RNorm ≤ ε
      inner_tired = restart ? inner_iter ≥ min(mem, inner_itmax) : inner_iter ≥ inner_itmax
      timer = time_ns() - start_time
      overtimed = timer > timemax_ns
      kdisplay(iter+inner_iter, verbose) && @printf("%5d  %5d  %7.1e  %.2fs\n", npass, iter+inner_iter, rNorm, ktimer(start_time))

      # Compute Vₖ₊₁.
      if !(solved || inner_tired || overtimed)
        if !restart && (inner_iter ≥ mem)
          push!(V, S(undef, n, p))
          push!(Z, S(undef, p, p))
        end
        V[inner_iter+1] .= Q
        Z[inner_iter+1] = Ztmp
      end
    end

    # Compute Yₖ by solving RₖYₖ = Zₖ with a backward substitution by block.
    Y = Z  # Yᵢ = Zᵢ
    for i = inner_iter : -1 : 1
      pos = nr + i - inner_iter         # position of Ψᵢ.ₖ
      for j = inner_iter : -1 : i+1
        mul!(Y[i], R[pos], Y[j], α, β)  # Yᵢ ← Yᵢ - ΨᵢⱼYⱼ
        pos = pos - j + 1               # position of Ψᵢ.ⱼ₋₁
      end
      Y[i] = Y[i] \ R[pos]  # Yᵢ ← Yᵢ \ Ψᵢᵢ
    end

    # Form Xₖ = NVₖYₖ
    for i = 1 : inner_iter
      mul!(Xr, V[i], Y[i], γ, β)
    end
    if !NisI
      solver.P .= Xr
      mulorldiv!(Xr, N, solver.P, ldiv)
    end
    restart && (X .+= Xr)

    # Update inner_itmax, iter, tired and overtimed variables.
    inner_itmax = inner_itmax - inner_iter
    iter = iter + inner_iter
    tired = iter ≥ itmax
    timer = time_ns() - start_time
    overtimed = timer > timemax_ns
  end
  (verbose > 0) && @printf("\n")

  # Termination status
  tired     && (status = "maximum number of iterations exceeded")
  solved    && (status = "solution good enough given atol and rtol")
  overtimed && (status = "time limit exceeded")

  # Update Xₖ
  warm_start && !restart && (X .+= ΔX)
  solver.warm_start = false

  # Update stats
  stats.niter = iter
  stats.solved = solved
  stats.timer = ktimer(start_time)
  stats.status = status
  return solver
end
