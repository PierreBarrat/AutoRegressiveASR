using DrWatson
@quickactivate("AutoRegressiveASR")

using AutoRegressiveASR
using Chain
using CSV
using DataFrames
using DataFramesMeta
using Dates
using DCATools
using DelimitedFiles
using JSON3
using StatsBase

function contact_prediction_data(basefolder::AbstractString)
    data = Dict{String,Any}("timestamp" => now())
    # Family specific: distances, contact prediction from equilibrium alignment, etc...
    ## will search the folder where the potts model is stored
    distances, scores_eq = get_distances_and_eq_scores(basefolder)
    data["distance_file"] = distances
    data["scores_plm_equilibrium_file"] = scores_eq
    data["ppv_eq"] = ppv(distances, scores_eq)

    # contact prediction from each tree
    distances = DataFrame(CSV.File(distances))
    data["asr"] = Dict()
    for repfol in ASRU.get_tree_folders(projectdir(basefolder, "data"))
        rep = @chain repfol splitpath last parse(Int, _)
        scores_ar, scores_iqtree, Meff_ardca, Meff_iqtree = get_scores_asr(repfol)
        data["asr"][rep] = Dict()
        data["asr"][rep]["ppv_iqtree"] = ppv(distances, scores_iqtree)
        data["asr"][rep]["Meff_iqtree"] = Meff_iqtree

        data["asr"][rep]["ppv_ardca"] = ppv(distances, scores_ar)
        data["asr"][rep]["Meff_ardca"] = Meff_ardca
    end
    @tag! data
    return data
end

function get_distances_and_eq_scores(basefolder)
    model_folder = @chain begin
        JSON3.read(joinpath(basefolder, "simulation_parameters.json"))
        _[:potts_file]
        dirname
        readdir(_; join=true)
    end
    distances = @chain begin
        findfirst(x -> occursin("_struct.csv", x), model_folder)
        @aside isnothing(_) && error("No struct file found model folder $(joinpath(basefolder, "simulation_parameters.json"))")
        model_folder[_]
    end
    scores_eq = @chain begin
        findfirst(x -> occursin("scores_plm_on_nat.csv", x), model_folder)
        @aside isnothing(_) && error("""
            No scores_plm_on_nat.csv file found model folder $(joinpath(basefolder, "simulation_parameters.json"))
            """)
        model_folder[_]
    end
    return distances, scores_eq
end

function get_scores_asr(folder)
    ar_scores = @chain begin
        joinpath(folder, "sample_root_plm/autoregressive/plm_scores.csv")
        CSV.File
        DataFrame
    end
    iqtree_scores = @chain begin
        joinpath(folder, "sample_root_plm/iqtree/plm_scores.csv")
        CSV.File
        DataFrame
    end

    Meff_ardca = @chain begin
        joinpath(folder, "sample_root_plm/autoregressive")
        X = readdir(; join=true)
        findfirst(x -> occursin(r"reconstructed_.*\.fasta", x), _)
        X[_]
        read_msa
        computeweights(; normalize=false)
        sum
    end
    Meff_iqtree = @chain begin
        joinpath(folder, "sample_root_plm/iqtree")
        X = readdir(; join=true)
        findfirst(x -> occursin(r"reconstructed_.*\.fasta", x), _)
        X[_]
        read_msa
        computeweights(; normalize=false)
        sum
    end

    return ar_scores, iqtree_scores, Meff_ardca, Meff_iqtree
end
