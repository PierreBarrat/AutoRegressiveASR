### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 80283bb4-d8c1-11ee-247b-d783b6ca3fbe
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")
	using AutoRegressiveASR
	using Chain
	using DataFrames
	using DataFramesMeta
	using DCATools
	using PlutoUI
	using StatsPlots
	using StatsBase
	using TreeTools
end

# ╔═╡ 7daed298-e497-4d88-bc3d-043dfe1f59bd
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ a10b7d75-df03-4a8f-9b11-82a11023f671
begin
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 4e1c78a5-80e6-4cee-90d9-fa3deb7f933d
PCF = pluto_ingredients(scriptsdir("figures_and_results/pca_functions.jl"))

# ╔═╡ dd73a005-3a35-46ad-82cb-119c538a1a84
strategies = PCF.strategies

# ╔═╡ 63c8922d-363d-47a8-bd6b-1b7696229c30
_fs = @bind folder_full Select(readdir(datadir("simulated/potts_yule"); join=true))

# ╔═╡ c9bdacf9-682b-461d-a467-9d45d69b82f8
path_folder, data_folder, prefix = let
	p, d = dirname(folder_full), basename(folder_full)
	prefix = split(d, "_")[1]
	p, d, prefix
end

# ╔═╡ 3c587c48-e5d6-4323-8508-ec477f5c2241


# ╔═╡ 2e57a7fe-0ac2-4556-9b35-7802bbd72bfe
data_all = let 
	config = Dict("folder" => folder_full)
	dat, _ = @produce_or_load(
		config;
		filename = x -> joinpath(x["folder"], "pca_data.jld2"), 
		suffix = "", 
		verbose=true, 
		tag=true, 
	) do config
		dat = PCF.pca_all_nodes(config["folder"], strategies)
	end
	dat
end

# ╔═╡ 2d6681a0-3155-400d-9d08-4f5f2eef6b19
function pca_nat_ref()
	p = plot()
	scatter!(
		data_all["equilibrium"][1,:], data_all["equilibrium"][2,:];
		marker = (3, .25, :lightblue, stroke(0)), label="",
	)
	plot!(
		xlabel = "PC-1",
		ylabel = "PC-2",
	)
	return p
end

# ╔═╡ df6f96ae-8147-49e6-bdb3-7abbda7dd716
begin
	pal = palette(:default)
	iqtree_color = pal[1]
	ar_color = pal[2]
	real_color = pal[3]

	real_marker = (32, :star)
	ml_marker = (18, :^)
	bayes_marker = (4, .75, stroke(0))
end

# ╔═╡ 37128352-bc1d-4c51-8057-b8d782d455f4
function strat_cosmetics(strat)
	color = if strat[1] == "iqtree"
		iqtree_color
	elseif strat[1] == "autoregressive"
		ar_color
	elseif strat[1] == "real"
		real_color
	end

	marker = if strat[1] == "real"
		real_marker
	elseif strat[2] == "ML"
		ml_marker
	elseif strat[2] == "Bayes"
		bayes_marker
	end

	label = if strat[1] == "real" || strat[2] == "ML"
		strat[1]
	else
		nothing
	end

	return color, marker, label
end

# ╔═╡ 4796d427-3716-4aa2-9fe3-4fd708cef2da
function plot_trace!(p, rep, list)
	for strat in (["real"], ["iqtree", "ML"], ["autoregressive", "ML"])
		x = map(n -> data_all[rep][n][strat][1,1], list)
		y = map(n -> data_all[rep][n][strat][2,1], list)
		color, _, _ = strat_cosmetics(strat)
		plot!(x, y; color, label="")
	end
	return p
end

# ╔═╡ a4a09a17-b1bc-4a1e-a1cb-9ac25b464241
function plot_strategy!(p, strat, X)
	color, marker, label = strat_cosmetics(strat)
	scatter!(p, [X[1,:]], [X[2,:]]; marker, color, label)
end

# ╔═╡ 4f40c0c4-922a-4933-8955-e3363cb160ba
rep = "16"

# ╔═╡ 907700af-44f5-4908-9c4c-2951ce2e8377
tree = read_tree(joinpath(folder_full, "data", rep, "tree.nwk"))

# ╔═╡ e5814bbd-21cb-46ab-bd9b-8bb44fffcca1
md"# Functions / Utils"

# ╔═╡ 0bd97b89-449a-447a-905f-944f929943c9
begin
	isbayes(strat) = (length(strat) > 1 && strat[2] == "Bayes")
	isml(strat) = (length(strat) > 1 && strat[2] == "ML")
	isreal(strat) = (strat[1] == "real")
end

# ╔═╡ cf78c311-5f4f-4a92-b2c0-9d53a24ac880
function plot_strategies!(
	p, rep, node;
	real = true,
	bayes = true,
	ml = true,
)
	@info strategies
	@info filter(s -> bayes || !isbayes(s), strategies)
	strats = @chain strategies begin
		filter(s -> real || !isreal(s), _)
		filter(s -> bayes || !isbayes(s), _)
		filter(s -> ml || !isml(s), _)
	end
	dat = data_all[rep][node]
	for strat in strats
		plot_strategy!(p, strat, dat[strat])
	end
	p
end

# ╔═╡ 56091e80-0dfd-49dc-98a0-3d7edaa1c729
function branch_with_most_nodes(tree)
	base = argmax(leaves(tree)) do leaf
		length(TreeTools.node_ancestor_list(leaf))
	end
	list = TreeTools.node_ancestor_list(base)[2:end-1]
	depths = map(node -> TreeTools.distance(tree, node, label(base)), list)
	return (nodes=list, depths=depths)
end

# ╔═╡ 4aa2b4d3-56db-4462-9102-b94a122dfa85
node_list = branch_with_most_nodes(tree).nodes

# ╔═╡ 7573d21a-dee0-4f44-83f4-9dde41f86fc4
let
	node = node_list[7]
	strat = strategies[1]
	p = pca_nat_ref()
	plot_strategies!(p, rep, node; bayes=false)
	# data_all[rep][node][strat]
	# strat
	plot_trace!(p, rep, node_list)
end

# ╔═╡ 81e1d628-7944-4a8f-8380-c0eea1c79129
let
	plts = [plot(rand(10)) for i in 1:10]
	@animate for p in plts
		p
	end
end

# ╔═╡ Cell order:
# ╠═80283bb4-d8c1-11ee-247b-d783b6ca3fbe
# ╠═7daed298-e497-4d88-bc3d-043dfe1f59bd
# ╠═a10b7d75-df03-4a8f-9b11-82a11023f671
# ╠═4e1c78a5-80e6-4cee-90d9-fa3deb7f933d
# ╠═dd73a005-3a35-46ad-82cb-119c538a1a84
# ╠═63c8922d-363d-47a8-bd6b-1b7696229c30
# ╠═c9bdacf9-682b-461d-a467-9d45d69b82f8
# ╠═3c587c48-e5d6-4323-8508-ec477f5c2241
# ╠═2e57a7fe-0ac2-4556-9b35-7802bbd72bfe
# ╠═2d6681a0-3155-400d-9d08-4f5f2eef6b19
# ╠═4796d427-3716-4aa2-9fe3-4fd708cef2da
# ╠═37128352-bc1d-4c51-8057-b8d782d455f4
# ╠═a4a09a17-b1bc-4a1e-a1cb-9ac25b464241
# ╠═cf78c311-5f4f-4a92-b2c0-9d53a24ac880
# ╠═df6f96ae-8147-49e6-bdb3-7abbda7dd716
# ╠═4f40c0c4-922a-4933-8955-e3363cb160ba
# ╠═7573d21a-dee0-4f44-83f4-9dde41f86fc4
# ╠═907700af-44f5-4908-9c4c-2951ce2e8377
# ╠═4aa2b4d3-56db-4462-9102-b94a122dfa85
# ╟─e5814bbd-21cb-46ab-bd9b-8bb44fffcca1
# ╠═0bd97b89-449a-447a-905f-944f929943c9
# ╠═56091e80-0dfd-49dc-98a0-3d7edaa1c729
# ╠═81e1d628-7944-4a8f-8380-c0eea1c79129
