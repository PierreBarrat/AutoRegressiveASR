using DrWatson
@quickactivate("AutoRegressiveASR")

using AncestralSequenceReconstruction
using ArDCA
using AutoRegressiveASR
using Dates
using JLD2
using JSON3

function asr_profile(parsed_args::AbstractDict; force=false)
    # Reconstruction strategies
    optbl = parsed_args["asr_opt_bl"]
    strategy_infer_bl = ASRMethod(; joint=false, optimize_branch_length_cycles=2, verbosity=2)
    strategy_ml = ASRMethod(; joint=false, ML=true, optimize_branch_length=false, verbosity=2)
    strategy_bayes = ASRMethod(;
        joint=false, ML=false, repetitions = 10, optimize_branch_length=false, verbosity=2,
    )
    prefix = "profile/"

    # Reconstruct on real tree using AR model
    @info "Reconstruction using profile"
    dat_folder = joinpath(parsed_args["folder"], "data")
    performed = false
    for fol in ASRU.get_tree_folders(dat_folder)
        @info fol

        evo_profile = ASR.ProfileModel(
            joinpath(fol, "alignment_leaves.fasta"), pc=0.1, reweighting=true,
        )

        if isdir(joinpath(fol, prefix))
            if force
                @warn "Removing $(joinpath(fol,prefix))"
                rm(joinpath(fol, prefix); recursive=true)
            else
                @warn "$(joinpath(fol, prefix)) already exists. Not running asr_profile again"
                continue
            end
        end
        mkpath(joinpath(fol, prefix));
        performed = true

        # reinfer branch length
        if Symbol(optbl) == :fromiqtree
            @info "Optimizing branch length starting from iqtree's tree"
            ASR.optimize_branch_length(
                joinpath(fol, "iqtree/tree_inferred.nwk"),
                joinpath(fol, "alignment_leaves.fasta"),
                evo_profile,
                strategy_infer_bl;
                outnewick = joinpath(fol, prefix, "tree_inferred.nwk"),
            )
        elseif Symbol(optbl) == :fromreal
            @info "Optimizing branch length starting from the real tree"
            ASR.optimize_branch_length(
                joinpath(fol, "tree.nwk"),
                joinpath(fol, "alignment_leaves.fasta"),
                evo_profile,
                strategy_infer_bl;
                outnewick = joinpath(fol, prefix, "tree_inferred.nwk"),
            )
        elseif Symbol(optbl) == :scale
            @info "Scaling branches of iqtree's tree"
            ASR.optimize_branch_scale(
                joinpath(fol, "iqtree/tree_inferred.nwk"),
                joinpath(fol, "alignment_leaves.fasta"),
                evo_profile,
                strategy_infer_bl;
                outnewick = joinpath(fol, prefix, "tree_inferred.nwk"),
            )
        elseif Symbol(optbl) == :real
            @info "Not optimizing branch length: using the real tree"
            cp(
                joinpath(fol, "tree.nwk"), joinpath(fol, prefix, "tree_inferred.nwk")
            )
        elseif Symbol(optbl) == :iqtree
            @info "Not optimizing branch length: using iqtree's tree"
            cp(
                joinpath(fol, "iqtree/tree_inferred.nwk"),
                joinpath(fol, prefix, "tree_inferred.nwk"),
            )
        else
            error("Unrecognized `optbl` $optbl")
        end

        # ML
        ASRU.reconstruct(
            fol, evo_profile, strategy_ml;
            tree_file = joinpath(prefix, "tree_inferred.nwk"),
            alignment_file = "alignment_leaves.fasta",
            outfiles = ["reconstructed_internals_ML.fasta"],
            prefix = prefix * "ML/"
        )
        # Bayesian
        ASRU.reconstruct(
            fol, evo_profile, strategy_bayes;
            tree_file = joinpath(prefix, "tree_inferred.nwk"),
            alignment_file = "alignment_leaves.fasta",
            outfiles = ["reconstructed_internals_rep$(i).fasta" for i in 1:strategy_bayes.repetitions],
            prefix = prefix * "Bayes/"
        )

    end

    if performed
        timestamp = now()
        parameters = @dict(
            arnet_file,
            timestamp,
            prefix,
            optbl,
            strategy_infer_bl,
            strategy_ml,
            strategy_bayes,
        )
        @tag!(parameters)
        open(joinpath(parsed_args["folder"], "profile_reconstruction_parameters.json"), "w") do f
            JSON3.pretty(f, JSON3.write(parameters))
        end
    end

    return dat_folder
end
