using Revise
using DrWatson
@quickactivate "AutoRegressiveASR"

using AutoRegressiveASR
using Chain
using Dates
using DelimitedFiles
using DCATools
using JLD2
using Random
using StatsBase

include(scriptsdir("families.jl"))

function simulate_dynamics(config)
    config = copy(config)
    id = savename(config["fam"]["prefix"], config; accesses = ["Ninit", "M"], sort = true)
    @info id
    folder = setup_folders!(config, config["fam"]; Ninit = config["Ninit"], id)
    # simulate for arnet
    tvals = simulate_arnet(config["fam"], folder, config["tvals_arnet"], config["M"])
    config["tvals_arnet_used"] = tvals

    # simulate for potts
    tvals = simulate_potts(config["fam"], folder, config["tvals_potts"], config["M"])
    config["tvals_potts_used"] = tvals

    config["timestamp"] = now()
    @tag!(config; force=true)
    save(joinpath(folder, "parameters.jld2"), config)
end

function setup_folders(config, fam; Ninit = 100, id = fam["prefix"])
    outfolder = datadir("dynamics", id)
    mkpath(outfolder)

    cp(fam["aln_nat"], joinpath(outfolder, "alignment_nat.fasta"); force=true)
    cp(fam["arnet_on_nat"], joinpath(outfolder, "arnet.jld2"); force=true)
    cp(fam["potts"], joinpath(outfolder, "potts.dat"); force=true)

    init_sequences = @chain "alignment_nat.fasta" begin
        joinpath(outfolder, _)
        X = read_msa
        randperm(length(_))[1:Ninit]
        DCATools.subsample(X, _)
    end
    write(joinpath(outfolder, "init_sequences.fasta"), init_sequences)

    config["aln_nat"] = joinpath(outfolder, "alignment_nat.fasta")
    config["arnet"] = joinpath(outfolder, "arnet.jld2")
    config["potts"] = joinpath(outfolder, "potts.dat")
    config["initial_sequences"] = joinpath(outfolder, "init_sequences.fasta")

    return outfolder
end

function simulate_arnet(fam, outfolder, tvals_in, M)
    @info "ArNet chains"
    arnet = JLD2.load(joinpath(outfolder, "arnet.jld2"))["arnet"]
    mkpath(joinpath(outfolder, "arnet_chains"))

    init_sequences = read_msa(joinpath(outfolder, "init_sequences.fasta"))

    tvals = unique(round.(tvals_in; sigdigits=3))
    if length(tvals) != length(tvals_in)
        @warn "Redundant time values (3 digits precision): $(collect(tvals_in)) vs $tvals"
    end
    writedlm(joinpath(outfolder, "arnet_chains/time_values.csv"), tvals)

    for (i, s0) in enumerate(init_sequences)
        print("$i ...")
        mkpath(joinpath(outfolder, "arnet_chains", "$i"))
        X_v_t = AutoRegressiveASR.propagate(s0, tvals, arnet, M)
        for (t, X) in zip(tvals, X_v_t)
            write(joinpath(outfolder, "arnet_chains/$(i)/alignment_t$(t).fasta"), X)
        end
    end
    println()

    return tvals
end

function simulate_potts(fam, outfolder, tvals_in, M)
    @info "Potts chain"
    potts = DCAGraph(joinpath(outfolder, "potts.dat"))
    mkpath(joinpath(outfolder, "potts_chains"))

    init_sequences = read_msa(joinpath(outfolder, "init_sequences.fasta"))

    all(isinteger, tvals_in) && @warn "time values should be sweeps, not swaps - got $tvals_in"
    tvals = @chain tvals_in*potts.L (@. round) (@. Int) unique

    # unique(round.(tvals_in; sigdigits=3))
    if length(tvals) != length(tvals_in)
        @warn """Redundant time values (3 digits precision): \
            Got: $(collect(tvals_in))
            When rounded: $tvals"""
    end
    writedlm(joinpath(outfolder, "potts_chains/time_values_swaps.csv"), tvals)

    for (i, s0) in enumerate(init_sequences)
        print("$i ...")
        mkpath(joinpath(outfolder, "potts_chains", "$i"))
        # swaps = Int.(round.(tvals * potts.L))
        X_v_t = AutoRegressiveASR.propagate(s0, tvals, potts, M) # Converting to mcmc swaps here
        for (t, X) in zip(tvals, X_v_t)
            write(joinpath(outfolder, "potts_chains/$(i)/alignment_t$(t).fasta"), X)
        end
    end
    println()

    return tvals
end

logrange(x, y, len; add_zero=false) = @chain begin
    range(log(x), log(y), length=len)
    @. exp
    collect
    add_zero ? vcat(0, _) : _
end

function tvals_potts(L)
    t1 = @chain range(0, L/4; length=10)/L collect
    t2 = logrange(5/16, 50, 41)
    return vcat(t1, t2)
end

function run()
    config = Dict(
       "tvals_arnet" => logrange(1e-2, 3, 51; add_zero=true),
       "tvals_potts" => logrange(1e-2, 40, 51; add_zero=true),
       "Ninit" => 100,
       "M" => 250,
    )

    for fam in values(families)
        p = copy(config)
        p["fam"] = fam
        simulate_dynamics(p)
    end
end

run()
