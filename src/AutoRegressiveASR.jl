module AutoRegressiveASR

using AncestralSequenceReconstruction
using ArDCA
using BackwardCoalescent
using Chain
using CSV
using Dates
using DataFrames
using DCATools
using DrWatson
using FASTX
using Smoothers
using StatsBase
using TreeTools

include("constants.jl")
export AA_IQTREE_ALPHABET

include("utils.jl")
export now_string, project_path

include("propagate_sequences.jl")

include("evolve_on_trees.jl")
export evolve, evolve_tree

include("simulation.jl")

include("ASRUtils/src/ASRUtils.jl")
const ASRU = ASRUtils
export ASRU, ASRUtils

end # module
