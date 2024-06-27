using DrWatson
quickactivate(@__DIR__, "AutoRegressiveASR")

using AncestralSequenceReconstruction
using BioSequenceMappings
using Dates
using JLD2
using JSON3
using Random
using StatsBase
using TreeTools

function pluto_ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
    m
end

STOOLS = pluto_ingredients(scriptsdir("DE_data/stiffler_data_functions.jl"))

const default_rec_params = Dict(
    "iqtree_model" => "JTT+R6",
    "arnet_file" => datadir("Stiffler/arnet/arnet_PF13354_lJ0.01_lH0.001.jld2"),
    "tree_style" => :star,
    "branch_opt" => :scale,
    "binarize" => true,
)

const default_data_params = Dict(
    "M" => 10,
    "tree_style" => :star,
    "seq_style" => :uniform,
    "suffix" => "",
)



#=RUN=#

data_params = copy(default_data_params)
rec_params = copy(default_rec_params)
do_cons = true
do_iqtree = true
do_arnet = true

#= Without binary tree (i.e. real star) =#
# data_params["suffix"] = "_nobinarize"
# rec_params["binarize"] = false
# do_cons = false
# do_iqtree = false

#= real star + opt branch length =#
# data_params["suffix"] = "_opt_nobin"
# rec_params["binarize"] = false
# rec_params["branch_opt"] = :opt
# do_cons = false
# do_iqtree = false

#= iqtree model finder =#
# data_params["suffix"] = "_modelfinder"
# rec_params["iqtree_model"] = ""
# do_cons = false
# do_arnet = false

#= arnet with gencode =#
data_params["suffix"] = "_gencode"
rec_params["with_code"] = true
do_cons = false
do_iqtree = false


# Mvals = 20 * (2 .^ (3:3))
# Mvals = [10, 40, 80]
Mvals = 20 * (2 .^ (3:5))
# Mvals = [160]
aln = read_fasta(datadir("Stiffler/aligned_data_ref/PSE1_rnd20_aligned_PF13354_noinserts.fasta"))
arnet = JLD2.load(default_rec_params["arnet_file"])["arnet"]
for M in Mvals, id in 1:100
    data_params["M"] = M
    folder, _ = STOOLS.select_data(data_params["M"], aln, data_params; id)
    STOOLS.reconstruct_all(folder, rec_params; arnet, do_iqtree, do_cons, do_arnet)
end

# for now this errors because I try to reroot iqtree's tree consistently with a star tree
# need to midpoint root probably
# # Tree: output of iqtree
# data_params["tree_style"] = :inferred
# rec_params["tree_style"] = :inferred
# folder, _ = STOOLS.select_data(data_params["M"], data_params)
# STOOLS.reconstruct_all(folder, rec_params)
