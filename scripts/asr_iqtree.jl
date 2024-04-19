using DrWatson
@quickactivate("AutoRegressiveASR")

using AutoRegressiveASR
using Dates
using JSON3

function asr_iqtree(parsed_args::AbstractDict; force=false)
    folder = parsed_args["folder"] |> abspath
    dat_folder = joinpath(folder, "data")

    # model = "Blosum62+I+G4" # model finder had this most of the time
    model = isempty(parsed_args["iqtree_model"]) ? nothing : parsed_args["iqtree_model"]
    bayes_repetitions = 5

    # Reconstruct using iqtree (will reinfer branches)
    performed = false
    @info "Reconstruction with iqtree"
    for fol in ASRU.get_tree_folders(dat_folder)
        dir_prefix = "iqtree/"
        iqtree_prefix = "IQTREE"
        if isdir(joinpath(fol, dir_prefix))
            if force
                @warn "Removing $(joinpath(fol,dir_prefix))"
                rm(joinpath(fol, dir_prefix); recursive=true)
            else
                @warn "$(joinpath(fol, dir_prefix)) already exists. Not running asr_iqtree again"
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
            joinpath(fol, dir_prefix); # operate in the fol/iqtree/ dir
            state_file = iqtree_prefix * ".state",
            prefix = "ML",
            out_files = ["reconstructed_internals_ML.fasta"],
            ML = true,
            alphabet = ASRU.AA_IQTREE_ALPHABET,
        )

        # Bayes - for testing
        ASRU.alignment_from_iqtree_state(
            joinpath(fol, dir_prefix); # operate in the fol/iqtree/ dir
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

    return dat_folder
end
