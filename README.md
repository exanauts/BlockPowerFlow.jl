# BlockPowerFlow

Prototype code to solve powerflow's linear system with multiple
right-hand side on the GPU.

Current, BlockPowerFlow implements:
- A block BICGSTAB method
- A wrapper to `CUSOLVERRF`

