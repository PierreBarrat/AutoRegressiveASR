### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ fd67489e-d881-11ee-3c76-c7e2a2b12431
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")
	using AutoRegressiveASR
	using Chain
	using CSV
	using DataFrames
	using DataFramesMeta
	using DCATools
	using MultivariateStats
	using StatsPlots
	using StatsBase
	using TreeTools
end

# ╔═╡ a84b764c-9e7f-4ce7-91ac-30810e314df9
base_folder = datadir(
	"simulated/potts_yule", 
	"PF00076_nleaves=50_nsweeps=15_ntrees=25_outgroup=false",
	"data/1/",
) 

# ╔═╡ 16d5484d-85aa-4016-9f99-a0801c93edd4
onehot = DCATools.onehot

# ╔═╡ d1ff144b-4a3f-4395-9dbe-0da68c004325
eq_sample = read_msa(joinpath(base_folder, "../../sample_adaBM_T200_PF00076.fasta"));

# ╔═╡ c40aa6ed-b318-4e13-b445-7f2ff5a335a7
reconstructed_nodes = @chain begin
	read_tree(joinpath(base_folder, "tree.nwk"))
	internals(; skiproot=true)
	map(label, _)
end

# ╔═╡ c4b41d34-6af2-4aa9-b4b0-f20a81ba4ebc
pca_model = fit(PCA, onehot(eq_sample); maxoutdim=2);

# ╔═╡ bd0611a3-19d1-4ec2-afe8-b61766935043
typeof(pca_model)

# ╔═╡ 44510d42-d2b2-4b88-b7ea-9073fb9c987a
    iqtree_profile = @chain begin
        joinpath(folder, "iqtree/IQTREE.state")
        ASRU.parse_iqtree_state_file
        getindex(_, target_node)
        getproperty(_, :model)
    end


# ╔═╡ 8f8da403-07b4-4b98-a7e7-12b9f1e4cf68
let
	strat = ["iqtree", "ML"]
	joinpath("here", strat..., "aagain")
end

# ╔═╡ d154e58d-e3fd-4a7e-9e2c-2e2a6b7ab1df
function pca_dat_folder(folder, nodes, pca_model::PCA)
	q = 21
	out = Dict()
	function pca_proj(strategy, node)
		is_iqtree = (splitpath(strategy)[1] == "iqtree")
		aln = if is_iqtree
			# profile = joinpath(folder, strategy, "IQTREE.state")
		else
			read_msa(joinpath(folder, strategy..., "reconstructed_$(node).fasta"))
		end
		return predict(pca_model, onehot(aln))
	end
	for node in nodes
		out[node] = Dict()
		# out[node]["pca_ar_bayes"] = pca_proj(["autoregressive", "Bayes"], node)
		# out[node]["pca_ar_ml"] = pca_proj("autoregressive_diversity/ML", node)
		# out["node"]["pca_iqtree_bayes"] = pca_proj("iqtree/Bayes", node)
		# out["node"]["pca_iqtree_ml"] = pca_proj("iqtree/Bayes", node)
	end
	return out
end

# ╔═╡ 60ddc2cd-6209-47dc-81e0-d0391445bb9c
pca_dat_folder(
	base_folder, reconstructed_nodes, pca_model
);

# ╔═╡ a0c54b34-1b27-4ef0-bfa3-0e95ef2c063e
md"# Utils"

# ╔═╡ f4f8f631-b5c9-4787-be2b-7d2cdc77c84b
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ e523f50b-fed4-4506-8159-79eed933a50f
DF = ingredients(scriptsdir("figures_and_results/diversity_functions.jl"))

# ╔═╡ 40e88159-4c7f-4d82-8d4a-aca4d193c94e
begin
	ml_from_profile = DF.ml_from_profile
	alignment_from_profile = DF.alignment_from_profile
end

# ╔═╡ 6f2d7e9f-70bd-420f-8b40-ab187f6c6a16
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

# ╔═╡ d6895bc4-48a2-4709-ae0a-a0848a93d0e2
S = alignment_from_strategy(base_folder, ["iqtree", "ML"], "internal_1")

# ╔═╡ Cell order:
# ╠═fd67489e-d881-11ee-3c76-c7e2a2b12431
# ╠═a84b764c-9e7f-4ce7-91ac-30810e314df9
# ╠═16d5484d-85aa-4016-9f99-a0801c93edd4
# ╠═d1ff144b-4a3f-4395-9dbe-0da68c004325
# ╠═c40aa6ed-b318-4e13-b445-7f2ff5a335a7
# ╠═c4b41d34-6af2-4aa9-b4b0-f20a81ba4ebc
# ╠═bd0611a3-19d1-4ec2-afe8-b61766935043
# ╠═60ddc2cd-6209-47dc-81e0-d0391445bb9c
# ╠═44510d42-d2b2-4b88-b7ea-9073fb9c987a
# ╠═8f8da403-07b4-4b98-a7e7-12b9f1e4cf68
# ╠═d154e58d-e3fd-4a7e-9e2c-2e2a6b7ab1df
# ╠═6f2d7e9f-70bd-420f-8b40-ab187f6c6a16
# ╠═d6895bc4-48a2-4709-ae0a-a0848a93d0e2
# ╟─a0c54b34-1b27-4ef0-bfa3-0e95ef2c063e
# ╠═e523f50b-fed4-4506-8159-79eed933a50f
# ╠═40e88159-4c7f-4d82-8d4a-aca4d193c94e
# ╠═f4f8f631-b5c9-4787-be2b-7d2cdc77c84b
