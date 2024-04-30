### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 0da80dce-fe43-11ee-22c0-43b4dadfc0b8
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

# ╔═╡ 8163f7c4-8b6e-42f9-b3ef-c368e800f8cb
begin
	local plot_defaults = joinpath(homedir(), ".julia/config/plot_defaults.jl")
	if isfile(plot_defaults)
		include(plot_defaults)
	end
end

# ╔═╡ 68bd14ba-77e6-4ca1-a3c6-7d8d15ce0221
md"# Setup"

# ╔═╡ 17c695dd-0fbc-413a-bda5-78b4afce93c8
if isdefined(@__MODULE__, :pubfig)
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 9888831f-2ccd-46ae-9f8a-0da660f1f36b
LM = pluto_ingredients(
	scriptsdir("figures_and_results/analyze_results_and_write_df.jl")
)

# ╔═╡ 6e82a4bc-2962-4262-a202-d372abc20d77
md"# Loading data"

# ╔═╡ 5c6b1a9b-000e-41d3-9463-0b8c22818148
basefolder = datadir("simulated/tree_regimes_arnet")

# ╔═╡ e7a25191-24f6-42c0-9ab1-259a3818f3b3
# ╠═╡ show_logs = false
data_all = let
	D = Dict()
	for folder in readdir(basefolder; join=true)
		@info folder
		params = parse_savename(basename(folder))[2]
		key = (params["coalescent"], params["treeheight"])
		@info key
		# Run measures on this simulation folder
		data, savefile = produce_or_load(
			Dict("folder" => folder);
			filename = x -> joinpath(x["folder"], "measures_asr"),
		) do config
			LM.analyze_results_and_write(config["folder"])
		end
		D[key] = data["asr"]
		@info "Saved measures in $savefile"
	end
	D
end;

# ╔═╡ f0c2d0e7-e2ac-4d8e-b594-08fdedf07a5c
md"# Figures"

# ╔═╡ b3c77e1c-57f6-43a7-9eb2-73025c0633c2
data_all |> values |> first |> values |> first |> names

# ╔═╡ 25e1c70d-5638-4c98-8266-92a848ed2504
md"## Branch length reconstruction"

# ╔═╡ 1ba0b6d9-946c-4d0f-b524-824a013d52e6
md"## Hamming distance to real"

# ╔═╡ 19abc312-0d49-4bd0-9b67-aa547d31e7a1
md"# Misc."

# ╔═╡ 6ce8902f-5630-4388-b335-55b6ea778b42
function sem(ystd, N)
	# error on mean calculated from sample standard deviation and number of samples
	# [-1.96, 196] has 95% of the mass for a normal distribution
	return 1.96 * ystd ./ sqrt.(N)
end

# ╔═╡ dcff23aa-bff6-4fbd-8906-dff729bf824e
begin
	# smoothing alg
	w = 20
	outliers_right = 0.
	smoothing_alg = :hist
end

# ╔═╡ 08ff2fd1-1b7b-40cc-b504-73cc02292144
function node_depth_series(data, strat)
	x, y, ystd, N = ASRU.easy_smooth(
		data[strat, "ML"], :node_depth, :node_depth_inferred; 
		w, alg=smoothing_alg, outliers_right
	)
	yerr = sem(ystd, N)

	return x, y, yerr
end

# ╔═╡ 44800b0f-35d0-49db-9fb0-7edd54107a55
function hamming_series(data, strat)
	x, y, ystd, N = ASRU.easy_smooth(
		data[strat, "ML"], :node_depth, :hamming_to_real_nogap; 
		w, alg=smoothing_alg, outliers_right
	)
	yerr = sem(ystd, N)

	return x, y, yerr
end

# ╔═╡ b29d094d-bc62-43ed-8fd6-5122d789ee16
regimes = data_all |> keys |> collect |> sort

# ╔═╡ ee63af55-11a8-4ac7-bef6-12823292f23d
begin
	yule(X) = filter(x -> x[1] ==("Yule"), X)
	kingman(X) = filter(x -> x[1] ==("Kingman"), X)
end

# ╔═╡ 3d0453ad-dc94-4e62-91a9-b4b60dad377d
regime_color = let
	pal = @chain regimes begin
		unique(x -> x[2], _)
		length
		palette(:bluesreds, _+2)
	end
	map(enumerate(regimes)) do (i, r)
		r => @chain length(regimes)/2 mod(i-1, _) Int pal[_+2]
	end |> Dict
end

# ╔═╡ 78ac8f04-6f0c-4727-9820-0b93296a6612
regime_line = Dict(
	s => (s[1] == "Kingman" ? :dash : :solid, regime_color[s])
	for s in regimes
) 

# ╔═╡ b2a2fbe8-048a-453f-a3e8-ac36888203e6
let p = plot()
	strat = "iqtree"
	for reg in yule(regimes)
		data = data_all[reg]
		x, y, yerr = node_depth_series(data, strat)
		plot!(
			p, x, y; 
			# ribbon = yerr, fillalpha=.2, fillcolor = regime_color[reg], 
			line = regime_line[reg], label=""
		)
	end
	plot!([0, 2], [0, 2]; line = (:black, :dashdot), label="")
	p
end

# ╔═╡ 1650dc4b-2652-43a3-a784-f6de4b84dd74
let p = plot()
	strat = "iqtree"
	for reg in yule(regimes)
		data = data_all[reg]
		x, y, yerr = hamming_series(data, strat)
		plot!(
			p, x, y; 
			# ribbon = yerr, fillalpha=.2, fillcolor = regime_color[reg], 
			line = regime_line[reg], label=""
		)
	end

	plot!(
		# xlims = (1e-3, 1.0),
		# xscale = :log10,
		ylim = (-0.025, 0.7),
	)
	
	p
end

# ╔═╡ Cell order:
# ╠═0da80dce-fe43-11ee-22c0-43b4dadfc0b8
# ╠═68bd14ba-77e6-4ca1-a3c6-7d8d15ce0221
# ╠═8163f7c4-8b6e-42f9-b3ef-c368e800f8cb
# ╠═17c695dd-0fbc-413a-bda5-78b4afce93c8
# ╠═9888831f-2ccd-46ae-9f8a-0da660f1f36b
# ╟─6e82a4bc-2962-4262-a202-d372abc20d77
# ╠═5c6b1a9b-000e-41d3-9463-0b8c22818148
# ╠═e7a25191-24f6-42c0-9ab1-259a3818f3b3
# ╟─f0c2d0e7-e2ac-4d8e-b594-08fdedf07a5c
# ╠═b3c77e1c-57f6-43a7-9eb2-73025c0633c2
# ╟─25e1c70d-5638-4c98-8266-92a848ed2504
# ╠═b2a2fbe8-048a-453f-a3e8-ac36888203e6
# ╠═08ff2fd1-1b7b-40cc-b504-73cc02292144
# ╟─1ba0b6d9-946c-4d0f-b524-824a013d52e6
# ╟─1650dc4b-2652-43a3-a784-f6de4b84dd74
# ╠═44800b0f-35d0-49db-9fb0-7edd54107a55
# ╟─19abc312-0d49-4bd0-9b67-aa547d31e7a1
# ╠═6ce8902f-5630-4388-b335-55b6ea778b42
# ╠═dcff23aa-bff6-4fbd-8906-dff729bf824e
# ╠═b29d094d-bc62-43ed-8fd6-5122d789ee16
# ╠═ee63af55-11a8-4ac7-bef6-12823292f23d
# ╠═3d0453ad-dc94-4e62-91a9-b4b60dad377d
# ╠═78ac8f04-6f0c-4727-9820-0b93296a6612
