# parsing command line arguments
include(joinpath(@__DIR__, "argparse.jl"))

const N = parsed_args["N"]
const HALFRADIUS = parsed_args["halfradius"]
const UPPER_BOUND = parsed_args["upper_bound"]

# Multithreading
# you may set this to the number of physical cores of your CPU
using LinearAlgebra
let nthr = min(4, Sys.CPU_THREADS ÷ 2)
    BLAS.set_num_threads(nthr)
    ENV["OMP_NUM_THREADS"] = nthr
end

# defining solvers/optimizers
include(joinpath(@__DIR__, "../src/optimizers.jl"))

include(joinpath(@__DIR__, "utils.jl"))

using Groups
import Groups.MatrixGroups

using PropertyT
import PropertyT.SA as StarAlgebras
import PropertyT.SW as SymbolicWedderburn
using PropertyT.PG # PermutationGroups

function wedderburn_decomposition(RG, Σ, action, psdrange)
    return SymbolicWedderburn.WedderburnDecomposition(
        Float64,
        Σ,
        action,
        StarAlgebras.basis(RG),
        StarAlgebras.Basis{UInt16}(@view StarAlgebras.basis(RG)[psdrange]),
    )
end