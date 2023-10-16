module BlockPowerFlow

using LinearAlgebra
using Printf
using SparseArrays

using KernelAbstractions

using CUDA
using CUDA.CUSPARSE
using CUDA.CUSOLVER

using AMDGPU
using AMDGPU.rocSPARSE
using AMDGPU.rocSOLVER

import LinearAlgebra.ldiv!

export CUSOLVERRF

include("preconditioners_nvidia.jl")
include("preconditioners_amd.jl")

include("block_utils.jl")
include("block_bicgstab.jl")
include("block_gmres.jl")

include("cusolverRF.jl")
using .CUSOLVERRF

end # module
