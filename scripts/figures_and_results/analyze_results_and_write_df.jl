using DrWatson
@quickactivate("AutoRegressiveASR")

using AutoRegressiveASR
using Chain
using CSV
using DataFrames
using DataFramesMeta
using Dates
using DCATools
using JLD2
using JSON3
using StatsBase

default_strategies = [
    ("iqtree", "ML"),
    ("iqtree", "Bayes"),
    ("autoregressive", "ML"),
    ("autoregressive", "Bayes"),
]

function analyze_results_and_write(folder; kwargs...)
    data = Dict{String, Any}("timestamp" => now())
    @tag!(data)

    # Measures of equilibrium properties of the teacher model
    parameters = JSON3.read(open(joinpath(folder, "simulation_parameters.json"), "r"))
    generative_model = load_generative_model(projectdir(parameters[:generative_model]))
    aln_eq = read_msa(parameters[:sample_equilibrium])

    data["likelihood_eq"] = map(s -> ASRU.myloglikelihood(s, generative_model), aln_eq)
    # data["pw_hamming_eq"] =

    # Measures of reconstruction quality for different strategies
    reconstruction_data = analyze_results(folder; kwargs...)
    data["asr"] = reconstruction_data

    return data
end



function analyze_results(
    folder::AbstractString;
    strategies = default_strategies,
)
    !isabspath(folder) && (folder = projectdir(folder))
    # retrieving generative model and sample
    parameters = JSON3.read(open(joinpath(folder, "simulation_parameters.json"), "r"))
    generative_model = load_generative_model(projectdir(parameters[:generative_model]))

    eq_sample_file = projectdir(parameters[:sample_equilibrium])
    sample_eq, model_consensus = if isfile(eq_sample_file)
        aln = read_msa(eq_sample_file)
        cons = DCATools.consensus(aln)
        cons = DCATools.num_to_aa(first(cons), cons.mapping)
        aln, cons
    else
        nothing, nothing
    end

    @info model_consensus

    data_folder = joinpath(folder, "data")

    # read data, measure observables and store results in dataframe
    data = Dict()
    for strat in strategies
        data[strat] = reconstruction_results(
            data_folder, strat; generative_model, model_consensus
        )
    end

    data[("real",)] = mapreduce(vcat, ASRU.get_tree_folders(data_folder)) do fol
        aln_leaves, aln_real, tree_real = real_files(fol)
        ASRU.reconstructed_to_df(
            tree_real, aln_leaves, aln_real;
            identifier = aln_real,
            alignment_real = aln_real,
            generative_model,
            model_consensus,
            inferred_tree_file = tree_real,
            include_reconstructed_sequence=true,
        )
    end

    return data
end

function real_files(folder)
    return (
        joinpath(folder, "alignment_leaves.fasta"),
        joinpath(folder, "alignment_internals.fasta"),
        joinpath(folder, "tree.nwk"),
    )
end
"""
    reconstructed_files(folder, strat)

`strat` must be like `("iqtree", "ML")`. Searches for trees in `folder/strat[1]` and for alignments in `folder/strat[1]/strat[2]`.
"""
function reconstructed_files(folder, strat)
    aln_files = @chain begin
        base = joinpath(folder, strat[1], strat[2])
        readdir
        filter(f -> occursin(r"\.fasta", f), _)
        joinpath.(base, _)
    end
    tree_file = joinpath(folder, strat[1], "tree_inferred.nwk")
    return aln_files, tree_file
end

function reconstruction_results(
    data_folder, strat;
    generative_model = nothing, model_consensus = nothing,
)
    # output is a DataFrame with all concatenating reconstructions for the given
    # strategy, accross folders in `data_folder`
    return mapreduce(vcat, ASRU.get_tree_folders(data_folder)) do fol
        aln_leaves, aln_real, tree_real = real_files(fol)
        alns, tree_inferred = reconstructed_files(fol, strat)
        # get all data for this folder
        df = mapreduce(vcat, alns) do aln
            ASRU.reconstructed_to_df(
                tree_real, aln_leaves, aln;
                identifier = aln,
                alignment_real = aln_real,
                generative_model,
                model_consensus,
                inferred_tree_file = tree_inferred,
                include_reconstructed_sequence=true,
            )
        end
        if length(groupby(df, [:label])) < size(df, 1)
            @combine groupby(df, [:label]) begin
                $(names(df, Number) .=> mean .=> identity)
            end
        else
            df
        end
    end
end




# function analyze_results_and_write_old(folder::AbstractString; kwargs...)
#     data_all = analyze_results(folder; kwargs...)

#     strats = collect(keys(data_all))
#     filenames = map(strats) do S
#         "data" * mapreduce(s -> "_"*s, *, S) * ".csv"
#     end

#     parameters = Dict(
#         "strategies" => strats,
#         "timesamp" => now(),
#         "filenames" => filenames,
#     )
#     @tag!(parameters)
#     JSON3.pretty(open(joinpath(folder, "analyze_results_log.json"), "w"), parameters)

#     for (fn, data) in zip(filenames, values(data_all))
#         CSV.write(joinpath(folder, fn), data)
#     end

#     return map(x -> joinpath(folder, x), filenames)
# end
