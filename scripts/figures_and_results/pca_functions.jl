using DrWatson
quickactivate(@__DIR__, "AutoRegressiveASR")
using AutoRegressiveASR
using Chain
using CSV
using DataFrames
using DataFramesMeta
using DCATools
using JSON3
using MultivariateStats
using TreeTools


strategies = [
    ["real"],
    ["iqtree", "Bayes"],
    ["autoregressive", "Bayes"],
    ["iqtree", "ML"],
    ["autoregressive", "ML"],
]

include(scriptsdir("figures_and_results/diversity_functions.jl"))

function pca_all_nodes(base_folder, strategies)
    # Common PCA model from an equilibrium sample of the Potts model
    sim_params = JSON3.read(joinpath(base_folder, "simulation_parameters.json"))
    eq_sample = read_msa(sim_params["sample_potts_file"])
    pca_model = fit(PCA, DCATools.onehot(eq_sample); maxoutdim=2);
    # for each sim in data/i, projection of all reconstructions
    data = Dict{String, Any}();
    proj_eq = predict(pca_model, DCATools.onehot(eq_sample))
    data["equilibrium"] = proj_eq
    for dat_folder in ASRU.get_tree_folders(joinpath(base_folder, "data"))
        id = splitpath(dat_folder)[end]
        data[id] = pca_dat_folder(dat_folder, strategies, pca_model)
    end
    return data
end

function pca_dat_folder(folder, strategies, pca_model)
    reconstructed_nodes = @chain begin
        read_tree(joinpath(folder, "tree.nwk"))
        internals(; skiproot=true)
        map(TreeTools.label, _)
    end
    return pca_dat_folder(folder, strategies, reconstructed_nodes, pca_model)
end

function pca_dat_folder(folder, strategies, nodes, pca_model::PCA)
    q = 21
    out = Dict()
    function pca_proj(strategy, node)
        aln = alignment_from_strategy(folder, strategy, node)
        return predict(pca_model, DCATools.onehot(aln))
    end
    for node in nodes
        out[node] = Dict()
        foreach(s -> out[node][s] = pca_proj(s, node), strategies)
    end
    leaves_real = read_msa(joinpath(folder, "alignment_leaves.fasta"))
    out["leaves"] = predict(pca_model, DCATools.onehot(leaves_real))
    return out
end

function alignment_from_strategy(folder, strategy, node)
    # folder is data/i/
    # strategy is ["iqtree", "ML"], or ["real"], ...
    q = 21
    return if strategy[1] == "real"
        aln = read_msa(joinpath(folder, "alignment_internals.fasta"))
        DCASample(aln[node])
    elseif strategy[1] == "autoregressive"
        read_msa(joinpath(
            folder,
            "autoregressive_diversity",
            strategy[2],
            "reconstructed_$(node).fasta"
        ))
    elseif strategy[1] == "iqtree"
        iqtree_state_file = joinpath(folder, "iqtree/IQTREE.state")
        profile = ASRU.parse_iqtree_state_file(iqtree_state_file)[node][:model]
        if strategy[2] == "ML"
            DCASample(ml_from_profile(profile))
        else
            alignment_from_profile(profile; M=500)
        end
    end
end
