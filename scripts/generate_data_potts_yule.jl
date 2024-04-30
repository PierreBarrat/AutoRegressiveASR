using DrWatson
@quickactivate "AutoRegressiveASR"

using AutoRegressiveASR
using BackwardCoalescent
using Dates
using DCATools
using JSON3
using TreeTools

function simulate_data_potts_yule(parsed_args::AbstractDict; force=false)
    ## Parameters
    # Potts model
    potts_file = parsed_args["generative_model"] |> abspath
    sample_potts_file = parsed_args["sample_equilibrium"]
    potts = DCAGraph(potts_file)

    # Genealogy
    nleaves = parsed_args["nleaves"] # number of leaves
    nsweeps = parsed_args["nsweeps"] # desired number of MCMC sweeps between root and leaves
    L = potts.L
    T = nsweeps * potts.L # desired height
    b = log(nleaves) / T # expected height log(nleaves)/b --> nsweeps * L
    outgroup = parsed_args["add_outgroup"]
    ntrees = parsed_args["ntrees"]
    nsim_per_tree = parsed_args["nsim_per_tree"]

    # Out folder
    timestamp = now()
    reps = nsim_per_tree
    generative_model = potts_file
    sample_equilibrium = sample_potts_file
    parameters = @dict(
        generative_model, sample_equilibrium,
        nleaves, nsweeps, ntrees, reps, b, outgroup,
        timestamp,
    )
    @tag!(parameters)
    identifier = savename(
        parsed_args["prefix"], parameters;
        accesses = [:nleaves, :nsweeps, :ntrees, :reps],
        sort = true,
    )
    outfolder = joinpath(parsed_args["outfolder"], identifier)

    if isdir(outfolder) || isfile(outfolder)
        if force
            @info "Remove all contents of $outfolder before simulating again, are you sure? [y/n]"
            yes = readline()
            !occursin("yes", yes) && (@warn "Aborting simulation"; return outfolder)
            try
                rm(outfolder; recursive=true)
            catch err
                @warn "Got error $err when trying to remove $outfolder"
            end
        else
            @warn "$outfolder already exists, not simulating again"
            return outfolder
        end
    end
    mkpath(outfolder)

    # Setting up folders
    @info "Using Potts model in $(project_path(potts_file))"
    @info "Results saved in $(project_path(outfolder))"

    if isfile(sample_potts_file)
        try
            cp(sample_potts_file, joinpath(outfolder, basename(sample_potts_file)))
        catch err
            @info """Tried to copy $sample_potts_file to $outfolder but got error $err
            Maybe the file already existed"""
        end
    elseif isempty(sample_potts_file)
        @info "No sample file given"
    else
        @warn "Could not find file $(abspath(sample_potts_file)) containing eq. sample of potts"
    end


    open(joinpath(outfolder, "simulation_parameters.json"), "w") do f
        JSON3.pretty(f, JSON3.write(parameters))
    end

    ## simulation
    dat_folder = joinpath(outfolder, "data")
    function get_tree()
        tree = genealogy(YuleCoalescent(nleaves, b))
        if parsed_args["normalize_tree_height"]
            Temp = TreeTools.distance_to_deepest_leaf(tree.root)
            foreach(nodes(tree; skiproot = true)) do n
                τ = branch_length(n)
                branch_length!(n, τ * T/Temp)
            end
        end
        return tree
    end

    generate_trees(
        dat_folder, get_tree;
        add_outgroup=outgroup,
        outgroup_distance = :auto,
        ntrees = ntrees,
        M = nsim_per_tree,
    )
    foreach(ASRU.get_tree_folders(dat_folder)) do fol
        AutoRegressiveASR.simulate_sequences(fol, potts)
    end

    @info "OUTPUT DIRECTORY $(abspath(outfolder))"
    return abspath(outfolder)
end

function generate_trees(
    outdir, rand_tree::Function;
    ntrees = 1,
    M = 1,
    add_outgroup=true,
    outgroup_distance = :auto,
    outgroup_name = "outgroup",
)
    cnt = 1
    for t in 1:ntrees
        tree_ref = rand_tree()
        for r in 1:M
            tree = copy(tree_ref)
            out = outdir * "/$(cnt)/"
            cnt += 1
            mkpath(out)
            # add outgroup
            if add_outgroup
                if outgroup_distance == :auto
                    outgroup_distance = mean(l -> distance(tree.root, l), leaves(tree))
                end
                graft!(tree, TreeNode(;tau = outgroup_distance, label=outgroup_name), tree.root)
            end
            TreeTools.label!(tree, tree.root, "root")
            # write
            write(out * "/tree.nwk", tree)
        end
    end
end

"""
    simulate_sequences(
        folder, model::ASR.EvolutionModel;
        tree_file, leaves_fasta, internals_fasta, prefix, kwargs...
    )

Use folder as base folder. Simulate sequences along `folder/tree_file` using `model`.
Store results as `folder/prefix/X_fasta` where `X` stands for leaves or internals.
Work is done by `ASRUtils.evolve`, additional `kwargs` are passed to it.
"""
function simulate_sequences(
    folder, model;
    tree_file = "tree.nwk",
    leaves_fasta = "alignment_leaves.fasta",
    internals_fasta = "alignment_internals.fasta",
    prefix = "",
    kwargs...
)
    tree = read_tree(folder * "/" * tree_file)
    return evolve(
        tree, model;
        leaves_fasta = folder * "/" * prefix *"/" * leaves_fasta,
        internals_fasta = folder * "/" * prefix *"/"* internals_fasta,
        kwargs...
    )
end
