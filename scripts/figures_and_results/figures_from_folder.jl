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

# ╔═╡ 157c8ab0-d63d-11ee-3794-27a66b691088
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

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

# ╔═╡ 7ebf97e2-fcca-4420-a01b-835de467ef89
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ c6d4a813-23f9-48a3-8418-398d94886da4
let
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 9faee986-bd06-477c-b0af-648394459cc4
plt_defaults = pubfig(20)

# ╔═╡ 367109e7-c4a5-46bf-8c9d-f4b6042c5cca
md"## Picking folder"

# ╔═╡ ecdbbfcf-0d17-49de-9f4f-18ae0d2e85eb
begin
	LM = pluto_ingredients(scriptsdir("figures_and_results/analyze_results_and_write_df.jl"))
	md"""
	This puts in scope the function `analyze_results_and_write(folder; strategies)`. It searches `folder` for `strat[1]/strat[2]` where `strat` is by default in 
	- `["iqtree", "ML"]
	- ...
	"""
end

# ╔═╡ 692ed38d-95b0-4003-b337-27c62efc9973
md"## Loading data"

# ╔═╡ c22d792b-fb09-4daa-9489-baab03481edf
cb_scatter = md"Scatter plots: $(@bind SP PlutoUI.CheckBox())"

# ╔═╡ e55196b3-a95e-4595-a601-e59115a477b5
cb_line = md"Smooth line plots: $(@bind LP PlutoUI.CheckBox())"

# ╔═╡ e0922d49-043e-494c-947e-e58e3a808224
md"# Figures"

# ╔═╡ e24faf1b-1034-4694-a8f3-9f7b173b867a
md"## Branch length"

# ╔═╡ 04ab3d51-0c30-4bd9-b08a-0c7e0540f3b1
md"## Hamming distance without gaps"

# ╔═╡ 0ab653c4-95e2-4676-bfda-d5b15f3d500f
_fs = @bind folder_full Select(readdir(datadir("simulated/potts_yule"); join=true))

# ╔═╡ f9b87a61-f1cd-4e57-a6c3-d2e6733a12ff
folder = basename(folder_full)

# ╔═╡ 6ebfeb60-837f-45e7-a918-bcbc91add5a2
title_id = split(folder, "_")[1]

# ╔═╡ f825cbb3-4ef9-4204-a0a4-d145be9cdd7b
data_all, _ = produce_or_load(
	Dict("folder" => folder_full);
	filename = x -> joinpath(x["folder"], "measures_asr.jld2"), suffix="",
) do config
	LM.analyze_results_and_write(config["folder"])
end;

# ╔═╡ b56973a7-2eb8-4805-a99a-6e180e25400b
data = data_all["asr"];

# ╔═╡ 39c9733b-6bbe-4d5d-bd70-3c90d5844b13
strategies = let
	st = collect(keys(data))
	
	lt(x,y) = if length(x) == length(y)
		x > y
	else
		length(x) > length(y)
	end
	sort(st; lt)
end

# ╔═╡ 577be925-516c-4efb-bbe0-3a07d6dfd2e0
simulation_parameters = JSON3.read(
	joinpath(folder_full, "simulation_parameters.json"), Dict
)

# ╔═╡ 39a98d0e-d550-4329-8048-b0462c8bd2e2
generative_model = simulation_parameters["potts_file"] |> DCAGraph

# ╔═╡ f4c24956-9f41-4b45-9617-ca965bcc7aab
model_consensus = let
	sample_file = projectdir(simulation_parameters["sample_potts_file"])
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

# ╔═╡ c71c1011-6c00-4947-ae82-2e0fa77d3446
figdir = mkpath(joinpath(plotsdir(), "potts_yule", basename(folder_full)))

# ╔═╡ 63b87289-1634-44ac-8e2c-e51ac3a7c271
md"## Hamming distance with gaps"

# ╔═╡ b96f65d5-cb63-4396-967a-f190e1979dcf
md"## Likelihood of reconstructed sequences"

# ╔═╡ 70b1c534-5f41-47f0-8fec-27201b0325d9
md"## Hamming distance to consensus"

# ╔═╡ 26f9e913-0b55-4d2f-9109-492787d386c6
md"## Functions & Misc"

# ╔═╡ 32ed2718-5cbc-4049-b66e-764f88aee027
begin
	# for histograms of likelihood
	deepnode = "internal_48"
end

# ╔═╡ 471f444a-7de5-4399-baba-8d189ed06670
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

# ╔═╡ 3fa3dd10-edbb-4ffe-94ec-55f665ba62b4
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

# ╔═╡ 085563cb-c3c4-4b28-bac6-ec58f79e4acf
let p = plot()
	for (i, strat) in enumerate(ml(strategies))
		x, y = ASRU.easy_smooth(
            data[strat], :node_depth, :node_depth_inferred; w, alg=smoothing_alg,
        )
		plot!(x, y, label=label2(strat), color=i, line=(3))
	end
	
	plot!(
		xlabel = "Real node depth (mcmc steps)",
		ylabel = "Inferred node depth",
		title = "$title_id - Inferred node depth",
		frame = :box,
		legend = :bottomright,
	)

	savefig(joinpath(figdir, "depth_inferred_v_true.png"))
	p
end

# ╔═╡ 8f6a4029-04e8-40f8-8dcd-144d5ede6c77
plt_hamming_to_real_nogaps = let p = plot()
	for strat in reconstruction(strategies)
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_real_nogap; w, alg=smoothing_alg,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real without gaps",
		frame = :box,
		legend = :bottomright,
	)

	savefig(joinpath(figdir, "hamming_to_real_ML_Bayes_nogaps.png"))
	p
end

# ╔═╡ c923e195-9628-4829-868b-dbb0ee6bf4dd
plt_hamming_to_real_nogaps_ml = let p = plot()
	for strat in ml(strategies)
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_real_nogap; w, alg=smoothing_alg,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real without gaps",
		frame = :box,
		legend = :bottomright,
	)

	savefig(joinpath(figdir, "hamming_to_real_ML_nogaps.png"))
	p
end

# ╔═╡ 6d13312e-f1cf-4d25-9eaa-6784ee6386f6
let p = plot()
	for (i, strat) in enumerate(bayesian(strategies))
		# @df data[strat] scatter!(
		# 	:node_depth, :hamming_to_real;
		# 	label=strat, color = i, marker=(3, .25, stroke(0))
		# )
		x, y = ASRU.easy_smooth(
            data[strat], :node_depth, :hamming_to_real; w, alg=smoothing_alg,
        )
		plot!(x, y, label=strat_label(strat), color=i, line=(3))
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real - Bayes",
		frame = :box,
		legend = :bottomright,
	)
end

# ╔═╡ 90467c05-4704-4903-8faa-076ab1461c4d
let p = plot()
	for (i, strat) in enumerate(ml(strategies))
		@df data[strat] scatter!(
			:node_depth, :hamming_to_real;
			label=label2(strat), color = i, marker=(6, .75, stroke(0))
		)
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real - ML",
		frame = :box,
		legend = :bottomright,
	)

	savefig(joinpath(figdir, "hamming_to_real_ML_scatter.png"))
	md"[Hamming to real ML scatter plot]"
end

# ╔═╡ 71613c45-8a7d-42ea-88bc-b4d5132f6a65
let p = plot()
	for (i, strat) in enumerate(ml(strategies))
		x, y = ASRU.easy_smooth(
            data[strat], :node_depth, :hamming_to_real; w, alg=smoothing_alg,
        )
		plot!(x, y, label=strat_label(strat), color=i, line=(3))
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real - ML",
		frame = :box,
		legend = :bottomright,
	)

	savefig(joinpath(figdir, "hamming_to_real_ML.png"))
	p
end

# ╔═╡ 103932ee-7c76-44f0-a5cb-78a141687996
plt_hamming_to_real_gaps = let p = plot()
	for strat in reconstruction(strategies)
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_real; w, alg=smoothing_alg,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real",
		frame = :box,
		legend = :bottomright,
	)

	savefig(joinpath(figdir, "hamming_to_real_ML_Bayes.png"))
	p
end

# ╔═╡ d13b857a-4071-4f96-816e-ea7949b3581c
plt_likelihood = let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :loglikelihood; w=15, alg=:hist,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Likelihood",
		title = "$title_id - Likelihood",
		frame = :box,
		legend = :topleft,
	)

	savefig(joinpath(figdir, "likelihood_v_depth.png"))

	p
end

# ╔═╡ c28c28dd-0e08-4e56-80d5-f0346b6ae6da
plt_likelihood_ml = let p = plot()
	for strat in vcat(ml(strategies), real(strategies))
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :loglikelihood; w=15, alg=:hist,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end

	hline!([mean(data_all["likelihood_eq"])], line = (2, :black, :dash), label="")
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Likelihood",
		title = "$title_id - Likelihood",
		frame = :box,
		legend = :topleft,
	)

	savefig(joinpath(figdir, "likelihood_v_depth_ML.png"))

	p
end

# ╔═╡ 41dbf0f5-8471-4d7d-901d-dd63b63783d2
plt_hamming_to_aln_consensus = let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_aln_consensus; w, alg=smoothing_alg,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
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

# ╔═╡ f0a26907-0304-4841-9401-9628a41a7c72
plt_hamming_to_model_consensus = let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_model_consensus; w, alg=smoothing_alg,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end

	# 
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "",
		title = "Hamming to generative model consensus",
		frame = :box,
		legend = :bottomleft,
	)

	savefig(joinpath(figdir, "hamming_to_model_consensus.png"))
	p
end

# ╔═╡ 2f154675-181a-40bf-9f2f-0b147cd4a4e2
reconstruction(strategies)

# ╔═╡ 7e189ca9-5bc9-402f-8f40-577a971d9016
real(strategies)

# ╔═╡ 1b4d2a5f-d121-47a2-a91f-f95a9dbd799e
function bar_plot_likelihood(data_all, wanted_strats)
	xlim = extrema(data_all["likelihood_eq"])

	# Equilibrium 
	p_dens = density(
		data_all["likelihood_eq"], line = (:black, :dash);
		label="", yticks = [], xlim, xticks = [],
		title = "Equilibrium",
		titlelocation = :left,
	)

	barplts = map(enumerate(wanted_strats)) do (i, strat)
		df = @subset data[strat] @byrow :label == deepnode
		p = bar(
			df.loglikelihood, ones(length(df.loglikelihood));
			label = "", 
			line = linestyle(strat), 
			color = strat_clr[strat],
			yticks = [], 
			xlim, 
			title = label_long(strat),
			titlelocation = :left,
		)
		if i < length(wanted_strats)
			plot!(p, xticks = [])
		else
			plot!(xlabel = "log-likelihood")
		end
		p
	end

	ν = .5
	L = length(barplts)
	return p = plot(
		p_dens, barplts...;
		layout = grid(1 + L, 1; heights=vcat(ν, ones(L)/L*(1-ν)))
	)
end

# ╔═╡ 4d76765d-db0d-465c-a56b-b75440dd2f0d
plt_bars_lk_ml = let
	wanted_strats = vcat(ml(strategies)..., ("real",))
	p = bar_plot_likelihood(data_all, wanted_strats)
	savefig(joinpath(figdir, "likelihood_at_root_ML.png"))
	p
end

# ╔═╡ 8ed5f31a-515d-4562-b42e-74673e2ae15a
plt_bars_lk_bayes = let
	wanted_strats = vcat(bayesian(strategies)..., ("real",))
	p = bar_plot_likelihood(data_all, wanted_strats)
	savefig(joinpath(figdir, "likelihood_at_root_Bayes.png"))
	p
end

# ╔═╡ Cell order:
# ╠═157c8ab0-d63d-11ee-3794-27a66b691088
# ╠═7ebf97e2-fcca-4420-a01b-835de467ef89
# ╠═c6d4a813-23f9-48a3-8418-398d94886da4
# ╠═9faee986-bd06-477c-b0af-648394459cc4
# ╟─367109e7-c4a5-46bf-8c9d-f4b6042c5cca
# ╠═f9b87a61-f1cd-4e57-a6c3-d2e6733a12ff
# ╠═6ebfeb60-837f-45e7-a918-bcbc91add5a2
# ╟─ecdbbfcf-0d17-49de-9f4f-18ae0d2e85eb
# ╟─692ed38d-95b0-4003-b337-27c62efc9973
# ╠═f825cbb3-4ef9-4204-a0a4-d145be9cdd7b
# ╠═b56973a7-2eb8-4805-a99a-6e180e25400b
# ╠═577be925-516c-4efb-bbe0-3a07d6dfd2e0
# ╠═39a98d0e-d550-4329-8048-b0462c8bd2e2
# ╠═f4c24956-9f41-4b45-9617-ca965bcc7aab
# ╠═39c9733b-6bbe-4d5d-bd70-3c90d5844b13
# ╠═c22d792b-fb09-4daa-9489-baab03481edf
# ╠═e55196b3-a95e-4595-a601-e59115a477b5
# ╟─e0922d49-043e-494c-947e-e58e3a808224
# ╠═c71c1011-6c00-4947-ae82-2e0fa77d3446
# ╟─e24faf1b-1034-4694-a8f3-9f7b173b867a
# ╟─085563cb-c3c4-4b28-bac6-ec58f79e4acf
# ╟─04ab3d51-0c30-4bd9-b08a-0c7e0540f3b1
# ╠═0ab653c4-95e2-4676-bfda-d5b15f3d500f
# ╟─8f6a4029-04e8-40f8-8dcd-144d5ede6c77
# ╟─c923e195-9628-4829-868b-dbb0ee6bf4dd
# ╟─63b87289-1634-44ac-8e2c-e51ac3a7c271
# ╟─6d13312e-f1cf-4d25-9eaa-6784ee6386f6
# ╟─90467c05-4704-4903-8faa-076ab1461c4d
# ╟─71613c45-8a7d-42ea-88bc-b4d5132f6a65
# ╟─103932ee-7c76-44f0-a5cb-78a141687996
# ╟─b96f65d5-cb63-4396-967a-f190e1979dcf
# ╟─d13b857a-4071-4f96-816e-ea7949b3581c
# ╟─c28c28dd-0e08-4e56-80d5-f0346b6ae6da
# ╟─4d76765d-db0d-465c-a56b-b75440dd2f0d
# ╟─8ed5f31a-515d-4562-b42e-74673e2ae15a
# ╟─70b1c534-5f41-47f0-8fec-27201b0325d9
# ╟─41dbf0f5-8471-4d7d-901d-dd63b63783d2
# ╟─f0a26907-0304-4841-9401-9628a41a7c72
# ╟─26f9e913-0b55-4d2f-9109-492787d386c6
# ╠═32ed2718-5cbc-4049-b66e-764f88aee027
# ╠═471f444a-7de5-4399-baba-8d189ed06670
# ╠═3fa3dd10-edbb-4ffe-94ec-55f665ba62b4
# ╠═2f154675-181a-40bf-9f2f-0b147cd4a4e2
# ╠═7e189ca9-5bc9-402f-8f40-577a971d9016
# ╠═1b4d2a5f-d121-47a2-a91f-f95a9dbd799e
