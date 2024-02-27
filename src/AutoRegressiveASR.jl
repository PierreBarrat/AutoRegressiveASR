module AutoRegressiveASR

using AncestralSequenceReconstruction
using ArDCA
using BackwardCoalescent
using Chain
using CSV
using Dates
using DataFrames
using DCATools
using FASTX
using Smoothers
using StatsBase
using TreeTools

include("constants.jl")
export AA_IQTREE_ALPHABET

include("utils.jl")
export now_string

include("propagate_sequences.jl")

include("evolve_on_trees.jl")
export evolve, evolve_tree

include("ASRUtils/src/ASRUtils.jl")
const ASRU = ASRUtils
export ASRU, ASRUtils

end # module
