### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 61726ef0-eb5a-11ee-1b76-4901ed3f59b2
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using AutoRegressiveASR
	using CSV
	using DataFrames
	using DataFramesMeta
	using DCATools
	using JLD2
	using JSON3
	using PlutoUI
	using StatsBase
	using StatsPlots
end

# ╔═╡ 0e693773-af4d-4155-a50e-ea8938a88d46
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 95de9723-9699-4168-84d7-ab4dc5844f5e
let
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ ac6a0e82-5dac-46c4-b871-dd1ece3d8d3d
plt_defaults = pubfig(14)

# ╔═╡ 30ad6744-ecfe-48d6-8628-252046ffc8a6
fam = "PF00076"

# ╔═╡ ff27ef68-53cb-4d22-87fb-a90b397f5436
split("PF00076_ad", r"-|_")

# ╔═╡ 152bf22f-5fe0-460c-901c-54a53bdb39d3
folders = @chain datadir("simulated/potts_yule") begin
	readdir(; join=true)
	filter(x -> split(basename(x), r"-|_")[1]==fam, _)
end

# ╔═╡ c5a2c7fd-9e17-430b-a9ea-ab4cd792e1ab
function key_from_folder(folder)
	params = joinpath(folder, "simulation_parameters.json")
	@assert isfile(params)

	params = JSON3.read(params, Dict)
	return (nsweeps = params["nsweeps"], nleaves = params["nleaves"])
end

# ╔═╡ 7ced3007-2a51-40df-8e12-ba09988ae7ec
data_all = let
	d = Dict()
	for folder in folders
		k = key_from_folder(folder)
		dat = JLD2.load(joinpath(folder, "measures_asr.jld2"))["asr"]
		d[k] = dat
	end
	d
end;

# ╔═╡ 1e61ba7a-444f-4516-9955-d6a4fb3eb63b
strategies = let
	st = collect(keys(first(values((data_all)))))
	
	lt(x,y) = if length(x) == length(y)
		x > y
	else
		length(x) > length(y)
	end
	sort(st; lt)
end

# ╔═╡ 0a922ce6-2be4-4201-954c-23761992dae4
param_groups = (
	cst_sweep = sort(filter(x -> x.nsweeps==10, collect(keys(data_all)))),
	cst_leaves = sort(filter(x -> x.nleaves==50, collect(keys(data_all)))),
)

# ╔═╡ 43434277-f088-4516-9842-89a667937e42
md"# Misc"

# ╔═╡ 793e5828-b6fe-4068-841b-b717d26faeba
begin
	# smoothing width
	w = 10
	smoothing_alg = :hist
	pal = palette(:default)
	strat_clr = Dict{Any,Any}(
		"iqtree" => pal[1], "autoregressive" => pal[2], "real" => pal[3]
	)
	for strat in strategies
		strat_clr[strat] = strat_clr[strat[1]]
	end
end

# ╔═╡ 80b262cf-066f-4dc7-a247-8226061fba31
let p = plot()
	strat = ("iqtree", "ML")
	for params in param_groups.cst_leaves
		data = data_all[params]
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_real_nogap; w, alg=smoothing_alg,
		)
		plot!(x, y, label="$params")
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$(strat) - Hamming to real without gaps",
		frame = :box,
		legend = :bottomright,
	)
	p
end

# ╔═╡ 3903bdd7-8af1-4c5d-9be5-8dd29de5a882
let p = plot()
	strat = ("autoregressive", "ML")
	for params in param_groups.cst_leaves
		data = data_all[params]
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_real_nogap; w, alg=smoothing_alg,
		)
		plot!(x, y, label="$params")
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$(strat) - Hamming to real without gaps",
		frame = :box,
		legend = :bottomright,
	)
	p
end

# ╔═╡ 7c7e8624-533a-4c9b-a51f-d47d83330b8c
let p = plot()
	strat_ar = ("autoregressive", "ML")
	strat_iqtree = ("iqtree", "ML")
	for params in param_groups.cst_leaves
		data = data_all[params]
		df_ar = @subset data[strat_ar] @byrow :label != "root"
		df_iq = data[strat_iqtree]
		x, y = ASRU.easy_smooth(
			df_ar.node_depth, # same for all strats, node depth in real tree
			df_iq.hamming_to_real_nogap .- df_ar.hamming_to_real_nogap;
			w, alg = smoothing_alg
		)

		plot!(x, y, label="$params")
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "ΔH to real",
		title = "Improvement: AR over iqtree",
		frame = :box,
		legend = :bottomleft,
	)
	p
end

# ╔═╡ f1e543d9-2023-4856-81b8-acd2ca9b1558
begin
	bayesian(strategies) = filter(x -> length(x)>1 && x[2]=="Bayes", strategies)
	ml(strategies) = filter(strategies) do x 
		length(x) < 2 && return false
		x[2] == "ML" || x[2] == "ml"
	end
	real(strategies) = filter(==(("real",)), strategies)
	reconstruction(strategies) = filter(!=(("real",)), strategies)
	strat_label(strat) = joinpath(strat...)

	iqtree(strategies) = filter(x -> x[1]=="iqtree", strategies)
	ar(strategies) = filter(x -> x[1]=="autoregressive", strategies)

	function label2(strat)
		length(strat) == 1 && return strat[1]
		strat[2] == "Bayes" ? "" : strat[1]
	end
	label_long(strat) = reduce((x,y) -> x*" - "*y, strat)

	
	function linestyle(strat)
		lw = 4
		return if length(strat) > 1 && strat[2] == "Bayes"
			(lw, :dash, strat_clr[strat[1]])
		else
			(lw, strat_clr[strat[1]])
		end
	end
	function barstyle(strat)
		(3, strat_clr[strat])
	end
end

# ╔═╡ Cell order:
# ╠═61726ef0-eb5a-11ee-1b76-4901ed3f59b2
# ╠═0e693773-af4d-4155-a50e-ea8938a88d46
# ╠═95de9723-9699-4168-84d7-ab4dc5844f5e
# ╠═ac6a0e82-5dac-46c4-b871-dd1ece3d8d3d
# ╠═30ad6744-ecfe-48d6-8628-252046ffc8a6
# ╠═ff27ef68-53cb-4d22-87fb-a90b397f5436
# ╠═152bf22f-5fe0-460c-901c-54a53bdb39d3
# ╠═c5a2c7fd-9e17-430b-a9ea-ab4cd792e1ab
# ╠═7ced3007-2a51-40df-8e12-ba09988ae7ec
# ╠═1e61ba7a-444f-4516-9955-d6a4fb3eb63b
# ╠═0a922ce6-2be4-4201-954c-23761992dae4
# ╟─80b262cf-066f-4dc7-a247-8226061fba31
# ╟─3903bdd7-8af1-4c5d-9be5-8dd29de5a882
# ╟─7c7e8624-533a-4c9b-a51f-d47d83330b8c
# ╠═43434277-f088-4516-9842-89a667937e42
# ╠═793e5828-b6fe-4068-841b-b717d26faeba
# ╠═f1e543d9-2023-4856-81b8-acd2ca9b1558
