using DrWatson
@quickactivate "AutoRegressiveASR"

using ArDCA
using AutoRegressiveASR
using BackwardCoalescent
using Dates
using Glob
using JLD2
using JSON3
using TreeTools

function simulate_data_arnet(parameters; force=false)
    @unpack ntrees, nleaves, treeheight, opt_bl, generative_model = parameters

    # Setting output folder
    identifier = savename(
        parameters["prefix"], parameters, parameters["suffix"];
        accesses = ["nleaves", "treeheight", "ntrees", "opt_bl"],
        sort = true,
    )
    outfolder = joinpath(parameters["outfolder"], identifier)
    ok = remove_outfolder_or_warn(outfolder, force)
    !ok && return outfolder
    mkpath(outfolder)
    @info "Results saved in $(projectdir(outfolder))"
    # saveing simulation parameters
    timestamp = now()
    log_parameters = @dict(
        generative_model,
        nleaves,
        treeheight,
        ntrees,
        opt_bl,
        timestamp,
    )
    @tag!(log_parameters)
    open(joinpath(outfolder, "simulation_parameters.json"), "w") do f
        JSON3.pretty(f, JSON3.write(log_parameters))
    end

    for rep in 1:ntrees
        dat_folder = joinpath(outfolder, "simulations/$(rep)/data")
        mkpath(dat_folder) # outfolder created automatically by this
        @info "Simulating tree $rep - results saved in $dat_folder"

        # Simulating trees
        @info "Sampling and re-rooting trees"
        simulate_trees(dat_folder, parameters)

        # Simulating sequences: only on the first tree
        @info "Using ArNet model in $(projectdir(generative_model))"
        @load parameters["generative_model"] arnet
        AutoRegressiveASR.simulate_sequences(
            joinpath(dat_folder), arnet;
            leaves_fasta = "alignment_leaves.fasta",
            internals_fasta = "alignment_internals.fasta",
        )
    end

    return abspath(outfolder)
end

function simulate_trees(outfolder, parameters)
    # Sample a tree that we will reroot n times
    tree = sample_random_tree(parameters)
    # The original tree (ultrametric) goes to outfolder
    # All re-rooted tree goes to outfolder/i (one of them identical to original)
    write(joinpath(outfolder, "tree.nwk"), tree)
    i = 1
    for n in traversal(tree, :preorder; leaves=false)
        mkpath(joinpath(outfolder, "$i"))
        root!(tree, label(n); remove_singletons=false)
        write(joinpath(outfolder, "$(i)", "tree.nwk"), tree)
        i += 1
    end
    return outfolder
end
function sample_random_tree(parameters)
    @unpack nleaves, treeheight = parameters
    tree = genealogy(YuleCoalescent(nleaves, 1))
    # Scaling branches
    treeheight_emp = TreeTools.distance_to_deepest_leaf(root(tree))
    foreach(nodes(tree; skiproot=true)) do n
        τ = branch_length(n)
        branch_length!(n, τ * treeheight/treeheight_emp)
    end

    return tree
end

function remove_outfolder_or_warn(outfolder, force)
    if isdir(outfolder) || isfile(outfolder)
        @info "$outfolder already exists"
        if force
            @info "Remove all contents of $outfolder before simulating again, are you sure? [yes/no]"
            yes = readline()
            if occursin("yes", yes)
                try
                    rm(outfolder; recursive=true)
                    return true
                catch err
                    @warn "Got error $err when trying to remove $outfolder"
                    return false
                end
            else
                @warn "Aborting simulation";
                return false
            end
        else
            @warn "$outfolder already exists, not simulating again"
            return false
        end
    end
    return true
end
