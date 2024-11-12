using DrWatson
@quickactivate "AutoRegressiveASR"

using AncestralSequenceReconstruction
using ArDCA
using Dates
using Glob
using JSON3
using JLD2

function asr_ardca(parameters::AbstractDict; force=false)
    @unpack dat_folder, generative_model, opt_bl = parameters

    if opt_bl != :real
        error("Expected `opt_bl==:real`, instead $(parameters["opt_bl"])")
    end

    # Preparing ardca model
    @load projectdir(generative_model) arnet
    evo_arnet = ASR.AutoRegressiveModel(arnet)

    # Reconstruction strategies
    strategy_ml = ASRMethod(; joint=false, ML=true, optimize_branch_length=false, verbosity=2)
    prefix = ""

    # Shared leaf alignment for all
    alignment_leaves = abspath(joinpath(dat_folder, "alignment_leaves.fasta"))

    # Reconstruct on real tree using AR model
    @info "Reconstruction using ArDCA"
    ## reconstruct on the original tree
    ASRU.reconstruct(
        dat_folder, evo_arnet, strategy_ml;
        tree_file = "tree.nwk",
        alignment_file = alignment_leaves, # is an abspath
        outfiles = ["original_reconstructed_internals_ML.fasta"],
        prefix,
    )
    performed = false
    for fol in glob([r"[0-9]+"], dat_folder)
        @info fol
        # check whether redo simulation or not
        if isfile(joinpath(fol, "reconstructed_internals_ML.fasta"))
            @warn "ArDCA reconstruction already performed for this file - skipping"
            continue
        end
        performed = true

        # ML
        ASRU.reconstruct(
            fol, evo_arnet, strategy_ml;
            tree_file = "tree.nwk",
            alignment_file = alignment_leaves, # is an abspath
            outfiles = ["reconstructed_internals_ML.fasta"],
            prefix,
        )
    end

    if performed
        timestamp = now()
        log_parameters = @strdict(
            timestamp,
            prefix,
            opt_bl,
            strategy_ml,
        )
        merge!(log_parameters, parameters)
        @tag!(log_parameters)
        open(joinpath(dat_folder, "ardca_reconstruction_parameters.json"), "w") do f
            JSON3.pretty(f, JSON3.write(parameters))
        end
    end

    return dat_folder
end
