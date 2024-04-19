using DrWatson
@quickactivate("AutoRegressiveASR")

using AncestralSequenceReconstruction
using ArDCA
using AutoRegressiveASR
using Dates
using IterTools
using JLD2
using JSON3

function asr_ardca_sample_internals(parsed_args::AbstractDict)
    return asr_ardca_sample_internals(;
        arnet_file = parsed_args["arnet"], folder = parsed_args["folder"]
    )
end

function asr_ardca_sample_internals(;
    arnet_file="", folder=nothing, force=false,
)
    dat_folder = joinpath(folder, "data")

    arnet = JLD2.load(projectdir(arnet_file))["arnet"]
    evo_arnet = ASR.AutoRegressiveModel(arnet)
    evo_profile = ASR.ProfileModel(arnet)

    # Reconstruction strategies
    opt_bl = :copy # the normal autoregressive run was done before, just get tree from there
    strategy_infer_bl = ASRMethod(; joint=false, optimize_branch_length_cycles=2, verbosity=2)
    strategy_ml = ASRMethod(; joint=false, ML=true, optimize_branch_length=false, verbosity=2)
    strategy_bayes = ASRMethod(;
        joint=false, ML=false, repetitions = 500, optimize_branch_length=false, verbosity=2,
    )
    prefix = "autoregressive_diversity/"

    timestamp = now()
    parameters = @dict(
        arnet_file,
        timestamp,
        prefix,
        opt_bl,
        strategy_infer_bl,
        strategy_ml,
        strategy_bayes,
    )
    @tag!(parameters)
    # will only be written if an actual computation is performed

    # Reconstruct on real tree using AR model
    @info "Reconstruction using ASR"
    performed = false
    for fol in takenth(ASRU.get_tree_folders(dat_folder), 5)
        @info fol
        if isdir(joinpath(fol, prefix))
            if force
                @warn "Removing $(joinpath(fol,prefix))"
                rm(joinpath(fol, prefix); recursive=true)
            else
                @warn "$(joinpath(fol, prefix)) already exists. Not running asr_ardca_sample_internals again"
                continue
            end
        end
        mkpath(joinpath(fol, prefix));

        # reinfer branch length
        if opt_bl == :opt
            @info "Optimizing branch length starting from iqtree's tree"
            ASR.optimize_branch_length(
                joinpath(fol, "iqtree/tree_inferred.nwk"),
                joinpath(fol, "alignment_leaves.fasta"),
                evo_profile,
                strategy_infer_bl;
                outnewick = joinpath(fol, prefix, "tree_inferred.nwk"),
            )
        elseif opt_bl == :scale
            @info "Optimizing branch scale using the original tree"
            ASR.optimize_branch_scale(
                joinpath(fol, "tree.nwk"),
                joinpath(fol, "alignment_leaves.fasta"),
                evo_profile,
                strategy_infer_bl;
                outnewick = joinpath(fol, prefix, "tree_inferred.nwk"),
            )
        elseif opt_bl == :copy
            cp(
                joinpath(fol, "autoregressive/tree_inferred.nwk"),
                joinpath(fol, prefix, "tree_inferred.nwk")
            )
        else
            @info "Not touching branch length"
        end

        # ML
        ASRU.reconstruct(
            fol, evo_arnet, strategy_ml;
            tree_file = joinpath(prefix, "tree_inferred.nwk"),
            alignment_file = "alignment_leaves.fasta",
            outfiles = "reconstructed.fasta",
            prefix = prefix * "ML/",
            alignment_per_node = true,
            node_list = nothing,
        )
        # Bayesian
        ASRU.reconstruct(
            fol, evo_arnet, strategy_bayes;
            tree_file = joinpath(prefix, "tree_inferred.nwk"),
            alignment_file = "alignment_leaves.fasta",
            outfiles = "reconstructed.fasta",
            outtable = "asr_table.csv",
            prefix = prefix * "Bayes/",
            alignment_per_node = true,
            node_list = nothing,
        )

        performed = true
    end

    if performed
        param_file = joinpath(folder, "ardca_diversity_reconstruction_parameters.json")
        open(param_file, "w") do f
            JSON3.pretty(f, JSON3.write(parameters))
        end
    end

    return dat_folder
end
