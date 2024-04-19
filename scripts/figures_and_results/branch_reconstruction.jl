using DrWatson
quickactivate(@__DIR__, "AutoRegressiveASR")

using AutoRegressiveASR
using Chain
using Combinatorics
using CSV
using DataFrames
using DataFramesMeta
using Dates
using IterTools
using JSON3
using StatsBase
using TreeTools

include(scriptsdir("families.jl"))

function measure_branch_reconstruction(basefolder::AbstractString)
    data = Dict{String,Any}("timestamp" => now())
    data_folder = joinpath(basefolder, "data")

    strategies = ["iqtree", "autoregressive"]

    # if potts used to simulate, measure real time in sweeps
    scale_factor = if occursin("potts", basefolder)
        fam = basefolder[findfirst(r"PF\d*", basefolder)]
        L = families[fam]["L"]
        @info "Scaling real branch length using length $L of family $fam"
        L
    else
        1
    end


    data_tmp = map(ASRU.get_tree_folders(data_folder)) do repfol
        node_data, pair_data = _measure_branch_reconstruction(repfol, strategies; scale_factor)
    end # array of tuple of `Dict("iqtree" => df, "autoregressive" => df)`

    data["nodes"] = Dict(
        strat => mapreduce(X -> X[1][strat], vcat, data_tmp) for strat in strategies
    )
    data["pairs"] = Dict(
        strat => mapreduce(X -> X[2][strat], vcat, data_tmp) for strat in strategies
    )
    @tag! data
    return data
end


function _measure_branch_reconstruction(folder::AbstractString, strategies; scale_factor=1)
    tree_real = read_tree(joinpath(folder, "tree.nwk"))
    tree_inf = Dict(s => read_tree(joinpath(folder, s, "tree_inferred.nwk")) for s in strategies)

    pair_data = measures_pairs(tree_real, tree_inf; scale_factor)
    branch_data = measures_branches(tree_real, tree_inf; scale_factor)
    return branch_data, pair_data
end

function measures_pairs(tree_real, tree_inf; scale_factor = 1)
    leaf_labels = map(TreeTools.label, leaves(tree_real))
    leaf_pairs = takenth(combinations(leaf_labels, 2), 10) |> collect

    dmat_real = map(leaf_pairs) do (x, y)
        distance(tree_real, x, y) / scale_factor
    end

    dat = Dict()
    for (strat, tree) in tree_inf
        dmat_inferred = map(x -> distance(tree, x[1], x[2]), leaf_pairs)
        dat[strat] = DataFrame(
            leaves = leaf_pairs,
            distance_real = dmat_real,
            distance_inferred = dmat_inferred,
        )
    end

    return dat
end

function measures_branches(tree_real, tree_inf; scale_factor = 1)
    node_labels = map(TreeTools.label, nodes(tree_real; skiproot=true))[1:2:end]
    bl_real = map(x -> branch_length(tree_real[x]) / scale_factor, node_labels)

    dat = Dict()
    for (strat, tree) in tree_inf
        bl_inf = map(x -> branch_length(tree[x]), node_labels)
        dat[strat] = DataFrame(
            node = node_labels,
            branch_length_real = bl_real,
            branch_length_inferred = bl_inf,
        )
    end

    return dat
end
