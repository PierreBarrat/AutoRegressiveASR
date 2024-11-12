using DrWatson
@quickactivate("AutoRegressiveASR")

using AutoRegressiveASR
using Dates
using JSON3

function asr_iqtree(parsed_args::AbstractDict; force=false)
    folder = parsed_args["folder"] |> abspath
    dat_folder = joinpath(folder, "data")


    model = isempty(parsed_args["iqtree_model"]) ? nothing : parsed_args["iqtree_model"]
    dir_prefix = dir_prefix_from_model(model)
    bayes_repetitions = 5

    # Reconstruct using iqtree (will reinfer branches)
    performed = false
    @info "Reconstruction with iqtree using model $model"
    for fol in ASRU.get_tree_folders(dat_folder)
        iqtree_folder = joinpath(fol, dir_prefix)
        iqtree_prefix = "IQTREE"
        if isdir(iqtree_folder)
            if force
                @warn "Removing $(iqtree_folder)"
                rm(iqtree_folder; recursive=true)
            else
                @warn "$(iqtree_folder) already exists. Not running asr_iqtree again"
                continue
            end
        end

        performed = true
        ASRU.reconstruct_iqtree(
            fol;
            tree_file = "tree.nwk",
            alignment_file = "alignment_leaves.fasta",
            out_tree_file = "tree_inferred.nwk", # ../ to store in folder and not in folder/prefix
            prefix = dir_prefix,
            iqtree_prefix,
            model,
        )

        # ML
        ASRU.alignment_from_iqtree_state(
            iqtree_folder; # operate in the fol/iqtree/ dir
            state_file = iqtree_prefix * ".state",
            prefix = "ML",
            out_files = ["reconstructed_internals_ML.fasta"],
            ML = true,
            alphabet = ASRU.AA_IQTREE_ALPHABET,
        )

        # Bayes - for testing
        ASRU.alignment_from_iqtree_state(
            iqtree_folder; # operate in the fol/iqtree/ dir
            state_file = iqtree_prefix * ".state",
            prefix = "Bayes",
            out_files = ["reconstructed_internals_rep$(i).fasta" for i in 1:bayes_repetitions],
            ML = false,
            alphabet = ASRU.AA_IQTREE_ALPHABET,
        )

        if get(parsed_args, "remove_iqtree_statefile", false)
            rm(joinpath(fol, dir_prefix, iqtree_prefix * ".state"))
        end
    end

    # Writing parameters
    if performed
        timestamp = now()
        parameters = @dict(
            model,
            bayes_repetitions,
            timestamp,
        )
        isnothing(parameters[:model]) && (parameters[:model] = :iqtree_model_finder)
        @tag!(parameters)
        open(joinpath(folder, "iqtree_reconstruction_parameters.json"), "w") do f
            JSON3.pretty(f, JSON3.write(parameters))
        end
    end

    return dat_folder, dir_prefix
end

dir_prefix_from_model(model) = isnothing(model) ? "iqtree/" : "iqtree-$(model)/"
