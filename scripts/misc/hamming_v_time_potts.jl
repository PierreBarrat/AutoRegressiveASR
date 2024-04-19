using Revise
using DrWatson
@quickactivate "AutoRegressiveASR"

using ArDCA
using AutoRegressiveASR
using Chain
using CSV
using DataFrames
using Dates
using DCATools
using JLD2
using StatsBase

include(scriptsdir("families.jl"))

# S ~ sample, s0 ~ init sequence
h_to_i(S, model, s0, t; kwargs...) = map(s -> DCATools.hamming(s, s0; normalize=true), S)
function loglk(S, model, s0, t; ref_model=model)
    return map(s -> round(ASRU.myloglikelihood(s, ref_model); sigdigits=2), S)
end
measures = (hamming_to_init = h_to_i, loglk = loglk)

function logrange(x, y, l; add_zero=false, round_to_int=true)
    R = @chain range(log(x), log(y); length=l) exp.(_)
    return collect(add_zero ? vcat(0, R) : R)
end
time_range_potts(L) = round.(Int, logrange(1, 250*L, 24; add_zero = true)) # in MCMC steps
time_range_exp = logrange(.01, 3, 24; add_zero=true)

function hamming_v_time_potts(fam; Nsamples = 100)
    potts = DCAGraph(fam["potts"])
    tvals = time_range_potts(potts.L)

    Seq = read_msa(fam["sample_potts_eq"])
    init_sequences = rand(Seq, Nsamples)
    # av_pw_distance_eq = mean(DCATools.pw_hamming_distance(Seq; step=100))

    vals = map(init_sequences) do s0
        AutoRegressiveASR.average_at_t(s0, measures, potts, tvals, 1)
    end

    df = DataFrame(t = tvals/potts.L, N = Nsamples * ones(Int, length(tvals)))
    for field in (:hamming_to_init, :loglk)
        dat = mapreduce(vcat, enumerate(tvals)) do (i, _)
            X = [getproperty(y[i], field) for y in vals]
            [mean(X) std(X)]
        end
        df[!, Symbol(field, :_av)] = dat[:, 1]
        df[!, Symbol(field, :_std)] = dat[:, 2]
        df[!, Symbol(field, :_conf95)] = 1.96 * df[!, Symbol(field, :_std)] / sqrt(Nsamples)
    end

    return df
end

function hamming_v_time_ardca(fam; Nsamples = 100)
    arnet = JLD2.load(fam["arnet"])["arnet"]
    tvals = time_range_exp

    Seq = read_msa(fam["sample_arnet_eq"])
    init_sequences = rand(Seq, Nsamples)
    # av_pw_distance_eq = mean(DCATools.pw_hamming_distance(Seq; step=100))

    vals = map(init_sequences) do s0
        AutoRegressiveASR.average_at_t(s0, measures, arnet, tvals, 1)
    end

    df = DataFrame(t = tvals, N = Nsamples * ones(Int, length(tvals)))
    for field in (:hamming_to_init, :loglk)
        dat = mapreduce(vcat, enumerate(tvals)) do (i, _)
            X = [getproperty(y[i], field) for y in vals]
            [mean(X) std(X)]
        end
        df[!, Symbol(field, :_av)] = dat[:, 1]
        df[!, Symbol(field, :_std)] = dat[:, 2]
        df[!, Symbol(field, :_conf95)] = 1.96 * df[!, Symbol(field, :_std)] / sqrt(Nsamples)
    end

    return df
end


savepath = datadir("dynamics")
for family in values(families)
    @info family
    produce_or_load(
        family, savepath;
        filename = X -> X["prefix"], prefix = "hamming_vs_time_potts_arnet", verbose=true
    ) do fam
        data_potts = hamming_v_time_potts(fam; Nsamples = 500)
        data_ardca = hamming_v_time_ardca(fam; Nsamples = 500)
        Dict(
            "family" => fam,
            "data_potts" => data_potts,
            "data_ardca" => data_ardca,
            "timestamp" => now()
        )
    end
end
