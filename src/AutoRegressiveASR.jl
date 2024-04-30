module AutoRegressiveASR

using AncestralSequenceReconstruction
using ArDCA
using BackwardCoalescent
import BioSequenceMappings: AbstractAlignment, pairwise_hamming
using Chain
using CSV
using Dates
using DataFrames
import DataFramesMeta: @orderby, @subset, @transform
using DCATools
using Distributions
using DrWatson
using FASTX
using JLD2
using OptimalTransport
using Tulip
using Smoothers
using StatsBase
using TreeTools

include("constants.jl")
export AA_IQTREE_ALPHABET

include("utils.jl")
export pluto_ingredients, project_path, load_generative_model

include("propagate_sequences.jl")

include("evolve_on_trees.jl")
export evolve, evolve_tree

include("simulation.jl")

include("dynamics_utils.jl")

include("ASRUtils/src/ASRUtils.jl")
const ASRU = ASRUtils
export ASRU, ASRUtils

include("contact_prediction.jl")
export ppv

include("distance.jl")

include("LBI.jl")
export lbi!

include("exp_distance_leaves.jl")
export weighted_leaf_distance!

end # module
