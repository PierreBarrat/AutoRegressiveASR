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

# ╔═╡ f9ae95d6-086f-11ef-3612-c55c7e9a9f18
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

# ╔═╡ 40134de3-c01b-4418-a229-4e077f5e39d4
begin
	local plot_defaults = joinpath(homedir(), ".julia/config/plot_defaults.jl")
	if isfile(plot_defaults)
		include(plot_defaults)
	end
end

# ╔═╡ 7b1c8ab6-89d8-4791-95a3-0c22dd6945b0
if isdefined(@__MODULE__, :pubfig)
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ e49257ff-afab-4bbc-8d01-55ed80820980
md"# Reading data"

# ╔═╡ e3f10a88-26ab-4e4d-8c81-d8c5d937cf64
md"## Folders"

# ╔═╡ d7080c4f-a482-4159-a272-965e57397641
folder_list = vcat(
	readdir(datadir("simulated/arnet_yule"); join=true),
)

# ╔═╡ 0c54ca8b-f3d6-4b74-b145-d009e84aa8cf
folder_picker = @bind folder_full Select(folder_list)

# ╔═╡ b2ac75ae-a46c-4b46-9e7b-1fe849208961
folder = basename(folder_full)

# ╔═╡ 24d6ab1c-9060-40ee-8745-a1bc38814142
figdir = mkpath(joinpath(
	plotsdir(), "paper/hamming/ardca_yule/", basename(folder_full)
))

# ╔═╡ 3461917c-e292-4db1-afe8-3a2052e325d9
md"## Reading data"

# ╔═╡ 246d3c2a-b706-43b2-a804-ef771d45e5d6
title_id = split(folder, "_")[1]

# ╔═╡ 398c97d1-eaa0-4187-b2e4-cd4873c3d7f2
HAM = pluto_ingredients(
	scriptsdir("figures_and_results/analyze_results_and_write_df.jl")
)

# ╔═╡ ec25bc9f-442e-45cf-8319-d550986d79a1
data_all, _ = produce_or_load(
	Dict("folder" => folder_full);
	filename = x -> joinpath(x["folder"], "measures_asr.jld2"), suffix="",
) do config
	HAM.analyze_results_and_write(config["folder"])
end;

# ╔═╡ e767d946-7364-4f3c-9855-3fae89baee32
data = data_all["asr"];

# ╔═╡ 436a6783-83f6-4ea0-b3c0-4312b0b2ed55
simulation_parameters = JSON3.read(
	joinpath(folder_full, "simulation_parameters.json"), Dict
)

# ╔═╡ 2dc59284-f5b3-4caf-a8b6-b9970431d0a9
generative_model = @chain begin
	simulation_parameters["generative_model"]
	load_generative_model
end

# ╔═╡ d814481b-a0e5-4f44-9500-e6247a874f40
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

# ╔═╡ 5d15c197-fdf0-447b-a662-96998303b362
strategies = let
	st = collect(keys(data))
	
	lt(x,y) = if length(x) == length(y)
		x > y
	else
		length(x) > length(y)
	end
	sort(st; lt)
end

# ╔═╡ ce1cc18f-e073-4157-86e0-b852cbdf474e
md"# Figures"

# ╔═╡ ec041156-eb28-4dc2-9d09-62d658bd4dc7
md"# Misc / Utils"

# ╔═╡ fad21c81-053a-4131-9fa6-f3286ecd18bf
function sem(ystd, N)
	# error on mean calculated from sample standard deviation and number of samples
	# [-1.96, 196] has 95% of the mass for a normal distribution
	return 1.96 * ystd ./ sqrt.(N)
end

# ╔═╡ 1eb4fd39-0b09-46e2-b2da-f7199746512c
begin
	# smoothing alg
	w = 20
	outliers_right = 0.
	smoothing_alg = :hist
end

# ╔═╡ 88a44666-cc56-4f43-addf-2dbc15ec58a5
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

# ╔═╡ 7e32bd49-0833-4a93-baa4-1ba7d9b33b7e
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

# ╔═╡ 75d9220d-46c4-4e11-9fc2-b1cea57b25bf
plt_likelihood_ml = let p = plot()
	for strat in vcat(ml(strategies), bayesian(strategies), real(strategies))
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :loglikelihood; 
			w=15, alg=:hist, outliers_right,
		)
		plot!(x, y, label=label_short(strat), line=linestyle(strat))
	end

	# hline!([mean(data_all["likelihood_eq"])], line = (2, :black, :dash), label="")

	lk_lim = let
		r, l = (1.15, 4)
		xmin, xmax = extrema(data_all["likelihood_eq"])
		mid = (xmax + xmin)/2
		mid - (mid - xmin)/l, mid + (xmax - mid)/r
	end
	density_lk_eq = let
		plt = density(data_all["likelihood_eq"]; xlim=lk_lim)
		x, y = plt[1][1][:x], plt[1][1][:y]
		plot(
			y, x; 
			ylim=lk_lim, label="", xticks = [], yticks = [], line=(3, :black),
			frame=:no, xaxis=false, fill=true, fillcolor = :black, fillalpha=.1
		)
	end
	
	plot!(p; 
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Likelihood",
		title = "$title_id - Likelihood",
		frame = :box,
		legend = :topleft,
		ylim = lk_lim,
	)

	savefig(joinpath(figdir, "likelihood_v_depth_ML.png"))

	p
	density_lk_eq

	plot(
		p, density_lk_eq;
		layout = @layout [a{0.91w} _ b{0.15w}]
	)
end

# ╔═╡ f7e153b6-270e-491b-80cb-c6a6561bdd0b
plt_hamming_to_aln_consensus = let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_aln_consensus; 
			w, alg=smoothing_alg, outliers_right,
		)
		plot!(x, y, label=label_short(strat), line=linestyle(strat))
	end

	# 
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "",
		title = "Hamming to alignment consensus",
		frame = :box,
		legend = :bottomleft,
	)

	savefig(joinpath(figdir, "hamming_to_aln_consensus.png"))
	p

end

# ╔═╡ Cell order:
# ╠═f9ae95d6-086f-11ef-3612-c55c7e9a9f18
# ╠═40134de3-c01b-4418-a229-4e077f5e39d4
# ╠═7b1c8ab6-89d8-4791-95a3-0c22dd6945b0
# ╠═e49257ff-afab-4bbc-8d01-55ed80820980
# ╠═e3f10a88-26ab-4e4d-8c81-d8c5d937cf64
# ╠═d7080c4f-a482-4159-a272-965e57397641
# ╠═b2ac75ae-a46c-4b46-9e7b-1fe849208961
# ╠═0c54ca8b-f3d6-4b74-b145-d009e84aa8cf
# ╠═24d6ab1c-9060-40ee-8745-a1bc38814142
# ╠═3461917c-e292-4db1-afe8-3a2052e325d9
# ╠═246d3c2a-b706-43b2-a804-ef771d45e5d6
# ╠═398c97d1-eaa0-4187-b2e4-cd4873c3d7f2
# ╠═ec25bc9f-442e-45cf-8319-d550986d79a1
# ╠═e767d946-7364-4f3c-9855-3fae89baee32
# ╠═436a6783-83f6-4ea0-b3c0-4312b0b2ed55
# ╠═2dc59284-f5b3-4caf-a8b6-b9970431d0a9
# ╠═d814481b-a0e5-4f44-9500-e6247a874f40
# ╠═5d15c197-fdf0-447b-a662-96998303b362
# ╠═ce1cc18f-e073-4157-86e0-b852cbdf474e
# ╟─75d9220d-46c4-4e11-9fc2-b1cea57b25bf
# ╟─f7e153b6-270e-491b-80cb-c6a6561bdd0b
# ╠═ec041156-eb28-4dc2-9d09-62d658bd4dc7
# ╠═fad21c81-053a-4131-9fa6-f3286ecd18bf
# ╠═1eb4fd39-0b09-46e2-b2da-f7199746512c
# ╠═88a44666-cc56-4f43-addf-2dbc15ec58a5
# ╠═7e32bd49-0833-4a93-baa4-1ba7d9b33b7e
