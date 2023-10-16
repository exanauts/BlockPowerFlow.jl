mutable struct AMD_IC0{SM}
  P::SM
end

function ic0(A::Union{ROCSparseMatrixCSC{T,Cint},ROCSparseMatrixCSR{T,Cint}}) where T <: LinearAlgebra.BlasFloat
  P = ic0(A, 'O')
  m,n = size(A)
  return AMD_IC0(P)
end

function ldiv!(y::ROCVector{T}, ic::AMD_IC0{ROCSparseMatrixCSR{T,Cint}}, x::ROCVector{T}) where T <: LinearAlgebra.BlasFloat
  copyto!(y, x)
  ldiv!(LowerTriangular(ic.P), y)   # Forward substitution with L
  ldiv!(LowerTriangular(ic.P)', y)  # Backward substitution with Lᴴ
  return y
end

function ldiv!(y::ROCMatrix{T}, ic::AMD_IC0{ROCSparseMatrixCSR{T,Cint}}, x::ROCMatrix{T}) where T <: LinearAlgebra.BlasFloat
  copyto!(y, x)
  ldiv!(LowerTriangular(ic.P), y)   # Forward substitution with L
  ldiv!(LowerTriangular(ic.P)', y)  # Backward substitution with Lᴴ
  return y
end

function ldiv!(y::ROCVector{T}, ic::AMD_IC0{ROCSparseMatrixCSC{T,Cint}}, x::ROCVector{T}) where T <: LinearAlgebra.BlasFloat
  copyto!(y, x)
  ldiv!(UpperTriangular(ic.P)', y)  # Forward substitution with L
  ldiv!(UpperTriangular(ic.P), y)   # Backward substitution with Lᴴ
  return y
end

function ldiv!(y::ROCMatrix{T}, ic::AMD_IC0{ROCSparseMatrixCSC{T,Cint}}, x::ROCMatrix{T}) where T <: LinearAlgebra.BlasFloat
  copyto!(y, x)
  ldiv!(UpperTriangular(ic.P)', y)  # Forward substitution with L
  ldiv!(UpperTriangular(ic.P), y)   # Backward substitution with Lᴴ
  return y
end

mutable struct AMD_ILU0{SM}
  P::SM
end

function ilu0(A::Union{ROCSparseMatrixCSC{T,Cint},ROCSparseMatrixCSR{T,Cint}}) where T <: LinearAlgebra.BlasFloat
  P = ilu0(A, 'O')
  return AMD_IC0(P)
end

function ldiv!(y::ROCVector{T}, ilu::AMD_ILU0{ROCSparseMatrixCSR{T,Cint}}, x::ROCVector{T}) where T <: LinearAlgebra.BlasFloat
  copyto!(y, x)
  ldiv!(UnitLowerTriangular(ilu.P), y)  # Forward substitution with L
  ldiv!(UpperTriangular(ilu.P), y)      # Backward substitution with U
  return y
end

function ldiv!(y::ROCMatrix{T}, ilu::AMD_ILU0{ROCSparseMatrixCSR{T,Cint}}, x::ROCMatrix{T}) where T <: LinearAlgebra.BlasFloat
  copyto!(y, x)
  ldiv!(UnitLowerTriangular(ilu.P), y)  # Forward substitution with L
  ldiv!(UpperTriangular(ilu.P), y)      # Backward substitution with U
  return y
end

function ldiv!(y::ROCVector{T}, ilu::AMD_ILU0{ROCSparseMatrixCSC{T,Cint}}, x::ROCVector{T}) where T <: LinearAlgebra.BlasFloat
  copyto!(y, x)
  ldiv!(LowerTriangular(ilu.P), y)      # Forward substitution with L
  ldiv!(UnitUpperTriangular(ilu.P), y)  # Backward substitution with U
  return y
end

function ldiv!(y::ROCMatrix{T}, ilu::AMD_ILU0{ROCSparseMatrixCSC{T,Cint}}, x::ROCMatrix{T}) where T <: LinearAlgebra.BlasFloat
  copyto!(y, x)
  ldiv!(LowerTriangular(ilu.P), y)      # Forward substitution with L
  ldiv!(UnitUpperTriangular(ilu.P), y)  # Backward substitution with U
  return y
end
