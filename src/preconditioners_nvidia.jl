mutable struct NVIDIA_IC0{SM,DM}
  P::SM
  z::DM
end

function ic0(A::Union{CuSparseMatrixCSC{T,Cint},CuSparseMatrixCSR{T,Cint}}, nrhs::Int) where T <: LinearAlgebra.BlasFloat
  P = ic02(A)
  m,n = size(A)
  z = nrhs == 1 ? CuVector{T}(undef, n) : CuMatrix{T}(undef, n, nrhs)
  return NVIDIA_IC0(P,z)
end

function ldiv!(y::CuVector{T}, ic::NVIDIA_IC0{<:CuVector{T},CuSparseMatrixCSR{T,Cint}}, x::CuVector{T}) where T <: LinearAlgebra.BlasFloat
  ldiv!(ic.z, LowerTriangular(ic.P), x)   # Forward substitution with L
  ldiv!(y, LowerTriangular(ic.P)', ic.z)  # Backward substitution with Lᴴ
  return y
end

function ldiv!(y::CuMatrix{T}, ic::NVIDIA_IC0{<:CuMatrix{T},CuSparseMatrixCSR{T,Cint}}, x::CuMatrix{T}) where T <: LinearAlgebra.BlasFloat
  ldiv!(ic.z, LowerTriangular(ic.P), x)   # Forward substitution with L
  ldiv!(y, LowerTriangular(ic.P)', ic.z)  # Backward substitution with Lᴴ
  return y
end

function ldiv!(y::CuVector{T}, ic::NVIDIA_IC0{<:CuVector{T},CuSparseMatrixCSC{T,Cint}}, x::CuVector{T}) where T <: LinearAlgebra.BlasFloat
  ldiv!(ic.z, UpperTriangular(ic.P)', x)  # Forward substitution with L
  ldiv!(y, UpperTriangular(ic.P), ic.z)   # Backward substitution with Lᴴ
  return y
end

function ldiv!(y::CuMatrix{T}, ic::NVIDIA_IC0{<:CuMatrix{T},CuSparseMatrixCSC{T,Cint}}, x::CuMatrix{T}) where T <: LinearAlgebra.BlasFloat
  ldiv!(ic.z, UpperTriangular(ic.P)', x)  # Forward substitution with L
  ldiv!(y, UpperTriangular(ic.P), ic.z)   # Backward substitution with Lᴴ
  return y
end

mutable struct NVIDIA_ILU0{SM,DM}
  P::SM
  z::DM
end

function ilu0(A::Union{CuSparseMatrixCSC{T,Cint},CuSparseMatrixCSR{T,Cint}}, nrhs::Int) where T <: LinearAlgebra.BlasFloat
  P = ilu02(A)
  m,n = size(A)
  z = nrhs == 1 ? CuVector{T}(undef, n) : CuMatrix{T}(undef, n, nrhs)
  return NVIDIA_ILU0(P,z)
end

function ldiv!(y::CuVector{T}, ilu::NVIDIA_ILU0{<:CuVector{T},CuSparseMatrixCSR{T,Cint}}, x::CuVector{T}) where T <: LinearAlgebra.BlasFloat
  ldiv!(ilu.z, UnitLowerTriangular(ilu.P), x)  # Forward substitution with L
  ldiv!(y, UpperTriangular(ilu.P), ilu.z)      # Backward substitution with U
  return y
end

function ldiv!(y::CuMatrix{T}, ilu::NVIDIA_ILU0{<:CuMatrix{T},CuSparseMatrixCSR{T,Cint}}, x::CuMatrix{T}) where T <: LinearAlgebra.BlasFloat
  ldiv!(ilu.z, UnitLowerTriangular(ilu.P), x)  # Forward substitution with L
  ldiv!(y, UpperTriangular(ilu.P), ilu.z)      # Backward substitution with U
  return y
end

function ldiv!(y::CuVector{T}, ilu::NVIDIA_ILU0{<:CuVector{T},CuSparseMatrixCSC{T,Cint}}, x::CuVector{T}) where T <: LinearAlgebra.BlasFloat
  ldiv!(ilu.z, LowerTriangular(ilu.P), x)      # Forward substitution with L
  ldiv!(y, UnitUpperTriangular(ilu.P), ilu.z)  # Backward substitution with U
  return y
end

function ldiv!(y::CuMatrix{T}, ilu::NVIDIA_ILU0{<:CuMatrix{T},CuSparseMatrixCSC{T,Cint}}, x::CuMatrix{T}) where T <: LinearAlgebra.BlasFloat
  ldiv!(ilu.z, LowerTriangular(ilu.P), x)      # Forward substitution with L
  ldiv!(y, UnitUpperTriangular(ilu.P), ilu.z)  # Backward substitution with U
  return y
end
