using DrWatson
quickactivate(@__DIR__, "AutoRegressiveASR")

using ArDCA
using AutoRegressiveASR
using BioSequenceMappings
using CSV
using DataFrames
using Dates
using JLD2
using JSON3
using Random
using StatsBase
using TreeTools

const non_gapped_positions = let
    wt_fasta = datadir("Stiffler/aligned_data_ref/PSE1_aligned_PF13354_noinserts.fasta")
    wt = read_fasta(wt_fasta) |> first |> collect
    findall(!=(1), wt)
end

function get_error_df(strat)
    base_folder = joinpath(datadir("Stiffler/subalignments/", strat))
    arnet = JLD2.load(datadir("Stiffler/arnet", "arnet_PF13354_lJ0.01_lH0.001.jld2"))["arnet"]

    Mvals = @chain readdir(base_folder) begin
       map(s -> match(r"M([0-9]*)", s), _)
       filter!(!isnothing, _)
       map(m -> parse(Int, m.captures[1]), _)
       sort
   end

   data = []
   for M in Mvals, folder in readdir(joinpath(base_folder, "M$M"); join=true)
       row = measure_row(folder; arnet)
       row["M"] = M
       row["id"] = splitpath(folder)[end]
       push!(data, row)
   end
   df = DataFrame(data)
   df = select(df, sort(names(df)))
   CSV.write(joinpath(base_folder, "data.csv"), df)
   return df
end

function measure_row(folder; arnet = nothing)
    seqs = get_sequences(folder)
    d = Dict{String, Any}()
    for label in (:cons, :iqtree, :arnet)
        #  Hamming distance to wild type
        Δpos, rec_aa, H = if ismissing(getfield(seqs, label))
            missing, missing, missing
        else
            pos = findall(non_gapped_positions) do i
                seqs.wt[i] != getfield(seqs, label)[i]
            end
            pos = map(i -> non_gapped_positions[i], pos)
            pos, [getfield(seqs, label)[i] for i in pos], length(pos)
        end

        d["pos_$label"] = Δpos
        d["rec_aa_$label"] = rec_aa
        d["nerr_$label"] = H

        # Likelihood in arnet
        d["likelihood_$(label)"] = if isnothing(arnet) || ismissing(getfield(seqs, label))
            missing
        else
            ArDCA.loglikelihood(getfield(seqs, label), arnet)
        end

    end
    return d
end

function get_sequences(folder)
    alphabet = Alphabet(:aa)
    wt_seq = wild_type(folder)
    cons_seq = rec_cons(folder)
    iqtree_seq = @chain rec_iqtree(folder) (ismissing(_) ? missing : alphabet(_))
    arnet_seq = rec_arnet(folder)
    return (wt = wt_seq, cons = cons_seq, arnet = arnet_seq, iqtree = iqtree_seq)
end

function rec_iqtree(folder)
    state_file = joinpath(folder, "iqtree/IQTREE.state")
    return if isfile(state_file)
        parsed_state_file = ASRU.parse_iqtree_state_file(state_file)
        parsed_state_file["root"].ml_seq
    else
        missing
    end
end
function rec_cons(folder)
    fasta = joinpath(folder, "consensus/consensus.fasta")
    return if isfile(fasta)
        @chain fasta read_fasta first collect
    else
        missing
    end
end
function rec_arnet(folder)
    files = readdir(folder; join=true)
    idx = findfirst(f -> occursin("arnet", f), files)
    isnothing(idx) && return missing

    fasta = joinpath(files[idx], "ML/ml.fasta")
    return if isfile(fasta)
        @chain fasta read_fasta find_sequence("root", _) collect(_[2])
    else
        missing
    end
end
wild_type(folder) = @chain joinpath(folder, "wt.fasta") read_fasta first collect
