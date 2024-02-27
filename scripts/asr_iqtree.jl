using DrWatson
@quickactivate("AutoRegressiveASR")

using AutoRegressiveASR
using JSON3

if isempty(ARGS) || length(ARGS) > 1
    println("Usage: `julia asr_iqtree.jl folder")
end

folder = ARGS[1] |> abspath
dat_folder = joinpath(folder, "data")

# Parameters
model = "Blosum62+I+G4" # model finder had this most of the time
# model = nothing # let's see what it finds
bayes_repetitions = 5

timestamp = now_string(; minute=true)
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


# Reconstruct using iqtree (will reinfer branches)
@info "Reconstruction with iqtree"
foreach(ASRU.get_tree_folders(dat_folder)) do fol
    dir_prefix = "iqtree/"
    iqtree_prefix = "IQTREE"

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
        out_files = ["reconstructed_internals_$(i).fasta" for i in 1:bayes_repetitions],
        ML = false,
        alphabet = ASRU.AA_IQTREE_ALPHABET,
    )
end
