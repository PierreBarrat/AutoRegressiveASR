using DrWatson
@quickactivate("AutoRegressiveASR")

using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "folder"
        help = "Folder containing data -- itself containing `data/i/tree.nwk` etc..."
        required = true
    "arnet"
        help = "JLD2 file with ardca model (in the `:arnet` field)"
        required = true
end
parsed_args = parse_args(ARGS, s)

using AncestralSequenceReconstruction
using ArDCA
using AutoRegressiveASR
using JLD2
using JSON3

arnet_file = project_path(parsed_args["arnet"])

arnet = JLD2.load(projectdir(arnet_file))["arnet"]
evo_arnet = ASR.AutoRegressiveModel(arnet)
evo_profile = ASR.ProfileModel(arnet)

# Reconstruction strategies
opt_bl = :opt
strategy_infer_bl = ASRMethod(; joint=false, optimize_branch_length_cycles=2, verbosity=2)
strategy_ml = ASRMethod(; joint=false, ML=true, optimize_branch_length=false, verbosity=2)
strategy_bayes = ASRMethod(;
    joint=false, ML=false, repetitions = 10, optimize_branch_length=false, verbosity=2,
)
prefix = "autoregressive/"

timestamp = now_string(; minute=true)
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
open(joinpath(parsed_args["folder"], "ardca_reconstruction_parameters.json"), "w") do f
    JSON3.pretty(f, JSON3.write(parameters))
end

# Reconstruct on real tree using AR model
@info "Reconstruction using ASR"
dat_folder = joinpath(parsed_args["folder"], "data")
foreach(ASRU.get_tree_folders(dat_folder)) do fol
    @info fol
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
    else
        @info "Not touching branch length"
    end

    # ML
    ASRU.reconstruct(
        fol, evo_arnet, strategy_ml;
        tree_file = joinpath(prefix, "tree_inferred.nwk"),
        alignment_file = "alignment_leaves.fasta",
        outfiles = ["reconstructed_internals_ML.fasta"],
        prefix = prefix * "ML/"
    )
    # Bayesian
    ASRU.reconstruct(
        fol, evo_arnet, strategy_bayes;
        tree_file = joinpath(prefix, "tree_inferred.nwk"),
        alignment_file = "alignment_leaves.fasta",
        outfiles = ["reconstructed_internals_$(i).fasta" for i in 1:strategy_bayes.repetitions],
        prefix = prefix * "Bayes/"
    )

end










