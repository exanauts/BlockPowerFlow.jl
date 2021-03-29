# BlockPowerFlow

Prototype code to solve powerflow's linear system with multiple
right-hand side on the GPU.

Currently, BlockPowerFlow implements:
- A block-BICGSTAB algorithm
- A wrapper to `CUSOLVERRF`

