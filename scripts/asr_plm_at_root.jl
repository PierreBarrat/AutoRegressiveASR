using DrWatson
@quickactivate("AutoRegressiveASR")

using AncestralSequenceReconstruction
using ArDCA
using AutoRegressiveASR
using Dates
using DelimitedFiles
using DCATools
using JLD2
using JSON3
using PlmDCA
using TreeTools

# to use alignment_from_profile for iqtree
include(scriptsdir("figures_and_results/diversity_functions.jl"))

function sample_internals_and_plm(args::AbstractDict)
    return sample_internals_and_plm(;
        arnet_file = args["arnet"], folder = args["folder"], iqtree_model = args["iqtree_model"],
    )
end

function get_target_node(folder::AbstractString)
    tree = read_tree(joinpath(folder, "tree.nwk"))
    # pick the child of the root closest to the root
    # not taking the root because iqtree does not reconstruct it (rerooting problems)
    n = argmin(branch_length, children(tree.root))
    return label(n)
end

function sample_internals_and_plm(;
    arnet_file="", iqtree_model="", folder=nothing, force=false,
)
    # results will be stored in data/i/prefix/
    # with subfolders for iqtree and autoregressive
    dat_folder = joinpath(folder, "data")
    prefix = "sample_root_plm/"


    # Number of samples at the internal node of interest
    M = 2_500
    weights = ones(M)/M

    # Reconstruction parameters for iqtree
    iqtree_model = isempty(iqtree_model) ? nothing : iqtree_model

    # Reconstruction parameters for ArDCA
    arnet = JLD2.load(projectdir(arnet_file))["arnet"]
    evo_arnet = ASR.AutoRegressiveModel(arnet)
    evo_profile = ASR.ProfileModel(arnet)
    strategy_infer_bl = ASRMethod(; joint=false, optimize_branch_length_cycles=2, verbosity=2)
    strategy_bayes = ASRMethod(;
        joint=false, ML=false, repetitions = M, optimize_branch_length=false, verbosity=2,
    )

    # parameters - written at the end if some simulation was done
    timestamp = now()
    parameters = @dict(
        iqtree_model,
        arnet_file,
        timestamp,
        prefix,
        strategy_infer_bl,
        strategy_bayes,
    )
    @tag!(parameters)


    # Reconstruction
    performed = false
    for fol in ASRU.get_tree_folders(dat_folder)
        @info fol
        outfolder = joinpath(fol, prefix)
        if isdir(outfolder)
            if force
                @warn "Removing $(outfolder) and simulating again"
                rm(outfolder; recursive=true)
            else
                @warn "$outfolder already exists. Delete folder manually to run again."
                continue
            end
        end
        mkpath(outfolder)

        # Label of target node (deepest of the tree that is not root)
        target_node = get_target_node(fol)

        # IQTREE
        iqtree_folder = joinpath(outfolder, "iqtree")
        iqtree_prefix = "IQTREE"
        ## this will generate  prefix/iqtree/IQTREE.state
        ASRU.reconstruct_iqtree(
            fol;
            tree_file = "tree.nwk",
            alignment_file = "alignment_leaves.fasta",
            out_tree_file = "tree_inferred.nwk", # ../ to store in folder and not in folder/prefix
            prefix = joinpath(prefix, "iqtree/"),
            iqtree_prefix,
            model = iqtree_model,
        )
        ## Read iqtree state file and sample for the right node
        state_file = joinpath(iqtree_folder, iqtree_prefix*".state")
        node_models = ASRU.parse_iqtree_state_file(
            state_file; alphabet=ASRU.AA_IQTREE_ALPHABET
        ) # Dict: label => (ml_seq=..., model=...)
        iqtree_aln = alignment_from_profile(node_models[target_node].model; M) # DCASample
        iqtree_out_file = joinpath(iqtree_folder, "reconstructed_$(target_node).fasta")
        @info iqtree_out_file
        write(iqtree_out_file, iqtree_aln)

        ## PLM on iqtree
        plm_iqtree = plmdca_asym(convert(Matrix, iqtree_aln.dat), weights)
        scores_iqtree = mapreduce(x -> collect(x)', vcat, plm_iqtree.score) # to have it as i j score matrix
        open(joinpath(iqtree_folder, "plm_scores.csv"), "w") do io
            writedlm(io, vcat(["i" "j" "score"], scores_iqtree), ' ')
        end

        # Ardca
        ardca_folder = joinpath(fol, prefix, "autoregressive")
        mkpath(ardca_folder)
        ASR.optimize_branch_length(
            joinpath(iqtree_folder, "tree_inferred.nwk"),
            joinpath(fol, "alignment_leaves.fasta"),
            evo_profile,
            strategy_infer_bl;
            outnewick = joinpath(ardca_folder, "tree_inferred.nwk"),
        )
        ASR.infer_ancestral(
            joinpath(ardca_folder, "tree_inferred.nwk"),
            joinpath(fol, "alignment_leaves.fasta"),
            evo_arnet,
            strategy_bayes;
            outfasta = joinpath(ardca_folder, "reconstructed.fasta"),
            alignment_per_node = true,
            node_list = [target_node],
        )
        ## plm ardca
        ardca_aln = read_msa(joinpath(ardca_folder, "reconstructed_$(target_node).fasta"))
        plm_ardca = plmdca_asym(convert(Matrix, ardca_aln.dat), weights)
        scores_ardca = mapreduce(x -> collect(x)', vcat, plm_ardca.score) # to have it as i j score matrix
        open(joinpath(ardca_folder, "plm_scores.csv"), "w") do io
            writedlm(io, vcat(["i" "j" "score"], scores_ardca), ' ')
        end

        #
        performed = true
    end

    param_file = joinpath(folder, "plm_at_root_log.json")
    open(param_file, "w") do f
        JSON3.pretty(f, JSON3.write(parameters))
    end

    return nothing
end
