### A Pluto.jl notebook ###
# v0.19.41

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

# ╔═╡ e3f4b0d4-f656-11ee-321f-d99da8673ea3
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using ArDCA
	using AutoRegressiveASR
	using CSV
	using DataFrames
	using DataFramesMeta
	using DCATools
	using JSON3
	using PlutoUI
	using StatsBase
	using StatsPlots
end

# ╔═╡ b6ad3d5f-c3a3-4e83-a1c4-b79ff68548bb
# ╠═╡ show_logs = false
begin
	local plot_defaults = joinpath(homedir(), ".julia/config/plot_defaults.jl")
	if isfile(plot_defaults)
		include(plot_defaults)
	end
end

# ╔═╡ 23867cb9-a0ab-46ab-8636-d58aa37ed5e8
PlutoUI.TableOfContents()

# ╔═╡ 1a359777-0e61-4bf7-bdd4-02c0e94fc77c
md"# Setup"

# ╔═╡ fea6a1f0-8ae5-4c4c-9f1c-ef7b74af3071
if isdefined(@__MODULE__, :pubfig)
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 87238856-61b3-42cb-ad07-e6248f79a3c5
LM = pluto_ingredients(
	scriptsdir("figures_and_results/analyze_results_and_write_df.jl")
)

# ╔═╡ bc4bb744-fe50-43bd-b4f9-857ce495cadd
md"# Reading data"

# ╔═╡ f78ea3ad-3686-4e75-9a98-174b5c6e07dd
md"## Folders"

# ╔═╡ d200f3c3-4b26-4345-b278-d10a7879eef8
folder_list = vcat(
	readdir(datadir("simulated/arnet_yule"); join=true),
)

# ╔═╡ e3aab824-e43a-425f-b151-eb32d0d07731
folder_picker = @bind folder_full Select(folder_list)

# ╔═╡ 807497b1-4ddd-471f-8c94-6d688c492c00
folder = basename(folder_full)

# ╔═╡ d0f90d56-f04c-456b-9ba7-27eb29262126
figdir = mkpath(joinpath(
	plotsdir(), "paper/hamming/ardca_yule/", basename(folder_full)
))

# ╔═╡ 17693eba-9a67-4483-90f6-7b0932726765
md"## Reading data"

# ╔═╡ 04dd0c96-4e78-4e9b-aa93-53aa6eb676b2
title_id = split(folder, "_")[1]

# ╔═╡ 715d6f56-365e-4d8f-bde8-4f02fd10d450
data_all, _ = produce_or_load(
	Dict("folder" => folder_full);
	filename = x -> joinpath(x["folder"], "measures_asr.jld2"), suffix="",
) do config
	LM.analyze_results_and_write(config["folder"])
end;

# ╔═╡ bcbf1882-221f-42f5-a004-b0cc888ace72
data = data_all["asr"];

# ╔═╡ 3e7bbd89-f2b8-4cc7-9474-011b2035e38a
simulation_parameters = JSON3.read(
	joinpath(folder_full, "simulation_parameters.json"), Dict
)

# ╔═╡ 78a77a70-3e2e-4234-b149-2d6506d42ae5
generative_model = @chain begin
	simulation_parameters["generative_model"]
	load_generative_model
end

# ╔═╡ 958ae778-ec6d-4930-b29c-c943cc3ded78
model_consensus = let
	sample_file = projectdir(simulation_parameters["sample_equilibrium"])
	Seq = if isnothing(sample_file) || !isfile(sample_file)
		T = 100
		M = 1000
		@info "sampling generative model for consensus. Check eq. time (default $T)"
		S = DCATools.sample(generative_model, M; Twait = T)
	else
		read_msa(sample_file)
	end
	cons = DCATools.consensus(Seq) # a DCASample object
	DCATools.num_to_aa(cons[1], cons.mapping) # a string sequence
end

# ╔═╡ 8e712428-e107-4f39-99d1-6b6fa2a08c66
strategies = let
	st = collect(keys(data))
	
	lt(x,y) = if length(x) == length(y)
		x > y
	else
		length(x) > length(y)
	end
	sort(st; lt)
end

# ╔═╡ 4e683b40-217f-4670-897a-0cfe6dc6185e
md"# Figures"

# ╔═╡ 6c8c0fc2-b6d1-4fbf-b5b4-5bf8a898a0cf
md"## Gaps - ML"

# ╔═╡ 53dca706-9b75-474a-8a1f-4010c6996dfe
folder_picker

# ╔═╡ 4648e611-d0de-4b23-b5f2-7b5463365999
md"# Misc / Utils"

# ╔═╡ 43e0f41e-b118-4d99-a753-2078e51afb86
function sem(ystd, N)
	# error on mean calculated from sample standard deviation and number of samples
	# [-1.96, 196] has 95% of the mass for a normal distribution
	return 1.96 * ystd ./ sqrt.(N)
end

# ╔═╡ 1de9ff83-6922-4d99-b723-8b45b8e97968
begin
	# smoothing alg
	w = 20
	outliers_right = 0.
	smoothing_alg = :hist
end

# ╔═╡ 7df463a0-e072-4ec2-82c3-350460c01837
begin
	# plot style
	pal = palette(:default)
	strat_clr = Dict{Any,Any}(
		"iqtree" => pal[1], "autoregressive" => pal[2], "real" => pal[3]
	)
	for strat in strategies
		strat_clr[strat] = strat_clr[strat[1]]
	end
end

# ╔═╡ 0daf2e97-bb2b-4e26-99a2-85e37071510e
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

	function label_short(strat)
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

# ╔═╡ 9de682ef-d1b5-4b43-9321-dbb946130df1
let p = plot()
	for (i, strat) in enumerate(ml(strategies))
		x, y, ystd, N = ASRU.easy_smooth(
            data[strat], :node_depth, :hamming_to_real; 
			w, alg=smoothing_alg, outliers_right
        )
		yerr = sem(ystd, N)
		plot!(
			x, y; ribbon = yerr, fillalpha=.2, 
			label=label_short(strat), line=linestyle(strat)
		)
	end

	# Difference
	S1, S2 = (("iqtree", "ML"), ("autoregressive", "ML"))
	D1 = sort(data[S1], :node_depth)
	D2 = sort(data[S2], :node_depth)
	X = D1.node_depth
	Y = D1.hamming_to_real - D2.hamming_to_real # iqtree - AR

	x, y, ystd, N = ASRU.easy_smooth(
		X, Y; w, alg=smoothing_alg, outliers_right, 
	)
	yerr = sem(ystd, N)
	plot!(
		x, y; ribbon = yerr, fillalpha=.2, label="improvement", color=:black
	)

	# 
	plot!(
		xlabel = "Node depth",
		ylabel = "Hamming distance to real",
		title = "",
		frame = :box,
		legend = :topleft,
	)

	savefig(joinpath(figdir, "hamming_to_real_gaps_ML.png"))
	p
end

# ╔═╡ bc99b3ef-d24a-48b3-8326-96baaab9fcbd
let p = plot()
	for (i, strat) in enumerate(ml(strategies))
		x, y, ystd, N = ASRU.easy_smooth(
            data[strat], :node_depth, :hamming_to_real_nogap; 
			w, alg=smoothing_alg, outliers_right
        )
		yerr = sem(ystd, N)
		plot!(
			x, y; ribbon = yerr, fillalpha=.2, 
			label=label_short(strat), line=linestyle(strat)
		)
	end

	# Difference
	S1, S2 = (("iqtree", "ML"), ("autoregressive", "ML"))
	D1 = sort(data[S1], :node_depth)
	D2 = sort(data[S2], :node_depth)
	X = D1.node_depth
	Y = D1.hamming_to_real_nogap - D2.hamming_to_real_nogap # iqtree - AR

	x, y, ystd, N = ASRU.easy_smooth(
		X, Y; w, alg=smoothing_alg, outliers_right, 
	)
	yerr = sem(ystd, N)
	plot!(
		x, y; ribbon = yerr, fillalpha=.2, label="improvement", color=:black,
	)

	# 
	plot!(
		xlabel = "Node depth",
		ylabel = "Hamming distance to real",
		title = "",
		frame = :box,
		legend = :topleft,
	)

	savefig(joinpath(figdir, "hamming_to_real_nogaps_ML.png"))
	p
end

# ╔═╡ 55be1196-a9df-40dc-992a-75fe0bc80d8e
let p = plot()
	for (i, strat) in enumerate(reconstruction(strategies))
		x, y, ystd, N = ASRU.easy_smooth(
            data[strat], :node_depth, :hamming_to_real_nogap; 
			w, alg=smoothing_alg, outliers_right
        )
		if strat[2] == "Bayes"
			yerr = sem(ystd, N)
			plot!(
				x, y; ribbon = yerr, fillalpha=.2, 
				label=strat[1], color = strat_clr[strat]
			)
		else
			plot!(
				x, y; label="", color = strat_clr[strat], line = (3, :dash)
			)
		end
	end

	# Difference
	S1, S2 = (("iqtree", "Bayes"), ("autoregressive", "Bayes"))
	D1 = sort(data[S1], :node_depth)
	D2 = sort(data[S2], :node_depth)
	X = D1.node_depth
	Y = D1.hamming_to_real_nogap - D2.hamming_to_real_nogap # iqtree - AR

	x, y, ystd, N = ASRU.easy_smooth(
		X, Y; w, alg=smoothing_alg, outliers_right, 
	)
	yerr = sem(ystd, N)
	plot!(
		x, y; ribbon = yerr, fillalpha=.2, label="improvement", color=:black
	)

	# Difference ML for ref
	S1, S2 = (("iqtree", "ML"), ("autoregressive", "ML"))
	D1 = sort(data[S1], :node_depth)
	D2 = sort(data[S2], :node_depth)
	X = D1.node_depth
	Y = D1.hamming_to_real_nogap - D2.hamming_to_real_nogap # iqtree - AR

	x, y, ystd, N = ASRU.easy_smooth(
		X, Y; w, alg=smoothing_alg, outliers_right, 
	)
	plot!(
		x, y; label="", line = (:black, :dash, 3)
	)

	# 
	plot!(
		xlabel = "Node depth",
		ylabel = "Hamming distance to real",
		title = "",
		frame = :box,
		legend = :topleft,
	)

	savefig(joinpath(figdir, "hamming_to_real_nogaps_Bayes.png"))
	p
end

# ╔═╡ Cell order:
# ╠═23867cb9-a0ab-46ab-8636-d58aa37ed5e8
# ╟─1a359777-0e61-4bf7-bdd4-02c0e94fc77c
# ╠═e3f4b0d4-f656-11ee-321f-d99da8673ea3
# ╠═b6ad3d5f-c3a3-4e83-a1c4-b79ff68548bb
# ╠═fea6a1f0-8ae5-4c4c-9f1c-ef7b74af3071
# ╠═87238856-61b3-42cb-ad07-e6248f79a3c5
# ╟─bc4bb744-fe50-43bd-b4f9-857ce495cadd
# ╟─f78ea3ad-3686-4e75-9a98-174b5c6e07dd
# ╠═d200f3c3-4b26-4345-b278-d10a7879eef8
# ╠═807497b1-4ddd-471f-8c94-6d688c492c00
# ╠═e3aab824-e43a-425f-b151-eb32d0d07731
# ╠═d0f90d56-f04c-456b-9ba7-27eb29262126
# ╟─17693eba-9a67-4483-90f6-7b0932726765
# ╠═04dd0c96-4e78-4e9b-aa93-53aa6eb676b2
# ╠═715d6f56-365e-4d8f-bde8-4f02fd10d450
# ╠═bcbf1882-221f-42f5-a004-b0cc888ace72
# ╠═3e7bbd89-f2b8-4cc7-9474-011b2035e38a
# ╠═78a77a70-3e2e-4234-b149-2d6506d42ae5
# ╠═958ae778-ec6d-4930-b29c-c943cc3ded78
# ╠═8e712428-e107-4f39-99d1-6b6fa2a08c66
# ╟─4e683b40-217f-4670-897a-0cfe6dc6185e
# ╟─6c8c0fc2-b6d1-4fbf-b5b4-5bf8a898a0cf
# ╠═9de682ef-d1b5-4b43-9321-dbb946130df1
# ╠═53dca706-9b75-474a-8a1f-4010c6996dfe
# ╟─bc99b3ef-d24a-48b3-8326-96baaab9fcbd
# ╟─55be1196-a9df-40dc-992a-75fe0bc80d8e
# ╟─4648e611-d0de-4b23-b5f2-7b5463365999
# ╠═43e0f41e-b118-4d99-a753-2078e51afb86
# ╠═1de9ff83-6922-4d99-b723-8b45b8e97968
# ╠═7df463a0-e072-4ec2-82c3-350460c01837
# ╠═0daf2e97-bb2b-4e26-99a2-85e37071510e
