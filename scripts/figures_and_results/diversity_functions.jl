using DrWatson
@quickactivate "AutoRegressiveASR"

using AutoRegressiveASR
using AncestralSequenceReconstruction
using Chain
using CSV
using DataFrames
using Dates
using DCATools
using JLD2
using StatsBase
using TreeTools

# Sample alignment from profile model
# Output a DCASample object, with the standard DCATools alphabet
# This is needed if the input profile model has a different alphabet (e.g. iqtree)
function alignment_from_profile(profile::ASR.ProfileModel; M = 1000)
    return @chain begin
        ASR.sample(profile, M)
        map(x -> ASR.intvec_to_sequence(x, profile.alphabet), eachcol(_))
        map(x -> ASR.sequence_to_intvec(x, ASR.Alphabet(:aa)), _)
        hcat(_...)'
        DCASample(_; mapping = ASR.reverse_mapping(ASR.Alphabet(:aa)))
    end
end

# ML sequence (Vector{Int}) from profile model
# Output a DCASample object, with the standard DCATools alphabet
# This is needed if the input profile model has a different alphabet (e.g. iqtree)
function ml_from_profile(profile::ASR.ProfileModel)
    return @chain begin ASR.ml_sequence(profile)
        ASR.intvec_to_sequence(profile.alphabet)
        ASR.sequence_to_intvec(ASR.Alphabet(:aa))
    end
end

# ╔═╡ 946c4e24-a011-4f60-9eec-de65b9042993
function _measures(profile::ASR.ProfileModel, ref_seqs::AbstractDict)
    reconstructed = alignment_from_profile(profile)
    M = distance_to_refseqs(reconstructed, ref_seqs)
    M["entropy"] = ASR.entropy(profile)
    return M
end
# ╔═╡ 72bf67cb-53b5-4595-8ea7-c3a1e34b8a55
function _measures(reconstructed::DCASample, rec_table, ref_seqs::AbstractDict)
    M = distance_to_refseqs(reconstructed, ref_seqs)
    M["entropy"] = -mean(rec_table.Total_LogLikelihood)
    return M
end

function distance_to_refseqs(reconstructed::DCASample, ref_seqs::AbstractDict)
    M = Dict{Any, Any}()
    # av and std hamming distance to set of ref. sequences (ml, real, consensus, ...)
    for (label, ref_seq) in ref_seqs
        m = distance_to_refseq(reconstructed, ref_seq, label)
        foreach(x -> M[x[1]] = x[2], m)
    end
    # self hamming distance
    self_hamming = DCATools.pw_hamming_distance(reconstructed; step=10)
    M["av_self_hamming"] = mean(self_hamming)
    M["std_self_hamming"] = std(self_hamming)

    return M
end

"""
"""
function distance_to_refseq(reconstructed::DCASample, ref_seq, label="real")
    hamming_to_ref = map(reconstructed) do s
        DCATools.hamming(s, ref_seq; normalize=true)
    end
    av_href, std_href = (mean(hamming_to_ref), std(hamming_to_ref))

    return Dict{Any,Any}(
        # "av_self_hamming" => av_self_hamming,
        # "std_self_hamming" => std_self_hamming,
        "av_hamming_$(label)" => av_href,
        "std_hamming_$(label)" => std_href,
    )
end

function diversity(folder::AbstractString, strat_folder, target_node::AbstractString, tree)

    leaves_real = read_msa(joinpath(folder, "alignment_leaves.fasta"))
    aln_consensus = DCATools.consensus(leaves_real)[1]

    internals_real = read_msa(joinpath(folder, "alignment_internals.fasta"))
    tarseq_real = internals_real[target_node]

    # ArDCA reconstruction (alignment + table with likelihood of reconstructed sequences)
    recseq_ardca = read_msa(
        joinpath(folder, strat_folder, "reconstructed_$(target_node).fasta")
    )
    mlseq_ardca = read_msa(joinpath(
        folder, split(strat_folder, '/')[1], "ML/reconstructed_$(target_node).fasta"
    ))[1]
    table_ardca = CSV.File(
        joinpath(folder, strat_folder, "asr_table_$(target_node).csv"),
    ) |> DataFrame

    # IQtree reconstruction: table with profile model at node of interest
    iqtree_profile = @chain begin
        joinpath(folder, "iqtree/IQTREE.state")
        ASRU.parse_iqtree_state_file
        getindex(_, target_node)
        getproperty(_, :model)
    end

    ref_seqs_iqtree = Dict(
        "ml" => ml_from_profile(iqtree_profile),
        "aln_consensus" => aln_consensus,
        "real" => tarseq_real,
    )
    ref_seqs_ardca = Dict(
        "ml" => mlseq_ardca,
        "aln_consensus" => aln_consensus,
        "real" => tarseq_real,
    )

    M_ardca = _measures(recseq_ardca, table_ardca, ref_seqs_ardca)
    M_iqtree = _measures(iqtree_profile, ref_seqs_iqtree)

    for M in (M_iqtree, M_ardca)
        M["name"] = target_node
        M["depth"] = TreeTools.distance_to_closest_leaf(tree, target_node)
        M["folder"] = folder
    end

    return DataFrame(M_ardca), DataFrame(M_iqtree)
end

function _diversity(folder::AbstractString, strategy = "autoregressive_diversity/Bayes")
    nodes = @chain begin
       readdir(joinpath(folder, strategy))
       filter(s -> occursin("reconstructed_", s), _)
       map(s -> match(r"reconstructed_(.*)\.fasta", s), _)
       map(m -> m.captures[1], _)
       filter(!=("root"), _)
   end

   tree = read_tree(joinpath(folder, "tree.nwk"))

   dat = map(node -> diversity(folder, strategy, node, tree), nodes)
   dat_asr = mapreduce(x -> x[1], vcat, dat)
   dat_iqtree = mapreduce(x -> x[2], vcat, dat)

   return (asr=dat_asr, iqtree=dat_iqtree)
end

function diversity_data(basefolder::AbstractString, outfile::AbstractString = "diversity_data.jld2")
    df_tuples = map(_diversity, ASRU.get_tree_folders(projectdir(basefolder, "data")))
    data = Dict{String,Any}("timestamp" => now())
    data["ardca"] = mapreduce(x -> x.asr, vcat, df_tuples)
    data["iqtree"] = mapreduce(x -> x.iqtree, vcat, df_tuples)
    @tag! data
    return data
end
