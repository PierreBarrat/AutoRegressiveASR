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

# ╔═╡ 81bc7c34-f1c5-11ee-17a1-f1983df9cbad
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using AutoRegressiveASR
	using CSV
	using DataFrames
	using DataFramesMeta
	using EasyFit
	using JSON3
	using Measures
	using PlutoUI
	using StatsBase
	using StatsPlots
end

# ╔═╡ ff4bc077-81bc-4dd7-80de-58ba09f4ff76
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ dea00e2e-8bea-4715-af4f-4b15ebeb792e
md"# Setup"

# ╔═╡ 4158b556-77a3-443c-96ef-2b90bab9ec97
let
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 93396504-09f3-44a5-a23d-834784b41218
plt_defaults = pubfig(20)

# ╔═╡ 4d06f815-17e9-41cf-9a5a-2c88e917665d
TOOLS = pluto_ingredients(scriptsdir(
	"figures_and_results/branch_reconstruction.jl"
))

# ╔═╡ 94e3fc55-8a73-494d-8e5e-8f7aae8b554d
fam_main = "PF00072"

# ╔═╡ 8ed7f4ca-343f-45d5-8333-25ee836a81db
savedir = projectdir("notes/article/figures/SI/")

# ╔═╡ 7a763db5-9890-4c8d-b468-c93239a2aeb4
fam_from_folder(folder) = split(basename(folder), "_")[1]

# ╔═╡ c9981050-3a9c-4cda-b1f3-145de001713e
method_from_folder(folder) = @chain folder begin
	splitpath
	_[end-1]
	split("_")
	getindex(1)
end

# ╔═╡ 6331fa1c-6d4c-43e2-94bd-d755f9a85e82
md"# Loading data"

# ╔═╡ 53b74403-c8b2-4bb0-9038-50ef587bc2bf
folder_list = vcat(
	readdir(datadir("simulated/potts_yule"); join=true),
	readdir(datadir("simulated/arnet_yule"); join=true),
)

# ╔═╡ 67bf7432-fd4d-410d-bb70-6e8344d7572e
_fs = @bind folder_full Select(folder_list)

# ╔═╡ ecd42c56-6100-4a99-a044-011bc93245db
begin
	folder = basename(folder_full)
	title_id = split(folder, "_")[1]
end

# ╔═╡ 6413baf8-deea-46e5-b609-896f95f7e4e8
branch_data, pair_data = let
	X, _ = produce_or_load(
		Dict("folder" => folder_full); 
		filename = x -> joinpath(x["folder"], "measures_branch_length"), 
	) do config
		TOOLS.measure_branch_reconstruction(config["folder"])
	end
	X["nodes"], X["pairs"]
end;

# ╔═╡ b93e7a60-9cbc-4f65-8083-ef58841ced47
simulation_parameters = JSON3.read(
	joinpath(folder_full, "simulation_parameters.json"), Dict
)

# ╔═╡ baee935b-9491-40d2-92eb-a648a190d1eb
strategies = collect(keys(branch_data))

# ╔═╡ 4a9f6f32-8ac5-43ee-8874-6a651546099f
max_dist_real = maximum(pair_data["iqtree"].distance_real)

# ╔═╡ 6ecb98ae-d792-4307-8a63-d212eb4c24fe
max_dist_inf = maximum(pair_data["autoregressive"].distance_inferred)

# ╔═╡ b587e52b-5330-4ced-9a40-2cb39dac0ae4
md"# Figures"

# ╔═╡ bd983086-a0ae-4072-b191-b1fcd63ed4e7
md"## Pairwise distance - Inferred vs real"

# ╔═╡ 1bc4405a-3249-4886-92e9-913b9516f221
names(pair_data["iqtree"])

# ╔═╡ fc4dabcc-dfe7-411b-a485-2559ab136b5a
_fs

# ╔═╡ 9389795a-45fb-4ea5-b43d-c0599a6fe218
md"## Pairwise distance - Cumulative distribution"

# ╔═╡ 6c407d72-2b03-4ae4-a95b-463fb438d711
md"## Error vs depth"

# ╔═╡ 439679e1-b121-4f82-a68e-a311af5280ee
md"# Misc"

# ╔═╡ 968c4424-fcc4-4d9e-b104-642e491e65f4
rel_err(a, ref) = (a-ref) /ref

# ╔═╡ 2484add9-eb75-4a32-81e6-6a5ac68a6b39
function bl_rel_error(data; threshold = 1e-1, max_err = Inf)
	return @chain data begin
		@subset(:branch_length_real .> threshold)
		@subset(:branch_length_inferred .> threshold)
		@transform(
			@byrow :rel_err_bl = rel_err(:branch_length_inferred, :branch_length_real)
		)
		@transform!(@byrow :rel_err_bl = min(max_err, :rel_err_bl))
	end
end

# ╔═╡ 1a2afd1a-6067-4a81-bda1-449c971d2060
begin
	# smoothing width
	w = 10
	smoothing_alg = :hist
end

# ╔═╡ a86493bd-daac-4302-b378-9c7c65db7152
# ╠═╡ disabled = true
#=╠═╡
begin
	iqtree(strategies) = filter(x -> x[1]=="iqtree", strategies)
	ar(strategies) = filter(x -> x[1]=="autoregressive", strategies)

	pal = palette(:default)
	strat_clr = Dict{Any,Any}(
		"iqtree" => pal[1], "autoregressive" => pal[2], "real" => pal[3]
	)
	for strat in strategies
		strat_clr[strat] = strat_clr[strat[1]]
	end
	
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
end
  ╠═╡ =#

# ╔═╡ 825d2389-d53c-4d2c-8c87-e521f050c5e0
begin
	local pal = palette(:default)
	strat_clr = Dict{Any,Any}(
		"iqtree" => pal[1], "autoregressive" => pal[2], "real" => pal[3]
	)

	function linestyle(strat)
		lw = 4
		return (lw, strat_clr[strat])
	end
end

# ╔═╡ 26699dfd-04f5-40ce-a59d-142603b5c1fe
plt_dist_inf_v_real = let p = plot()
	for strat in strategies
		@df pair_data[strat] scatter!(
			:distance_real, :distance_inferred;
			label = "", marker = (3, .5, stroke(0), strat_clr[strat]),
		)
	end
	
	# Linear fit AR v real using the f leftmost points
	f = 0.05
	X, Y = @chain pair_data["autoregressive"] begin
		sort(:distance_real)
		_.distance_real, _.distance_inferred
	end

	if method_from_folder(folder_full) == "arnet"
		plot!([0, max_dist_real], [0, max_dist_real], label="", line = (:black, :dash))
	end
	
	# Linear fit on short distances
	# L = Int(round(length(X) * f))
	# linfit = fitlinear(X[1:L], Y[1:L])
	# plot!(X, linfit.(X), line=(:black), label="")
	
	plot!(
		xlabel = "real distance",
		ylabel = "inferred distance",
		title = "Distance between leaves",
		ylims = (-0.25, max_dist_inf * 0.7),
	)
end

# ╔═╡ ab30475e-a052-4dc4-b316-42b33d2d2205
plt_dist_cumulative = let p = plot()
	dvals = range(0, max_dist_real * 1.2, length=1000)

	
	for strat in strategies
		cdf = ecdf(pair_data[strat].distance_inferred)
		plot!(dvals, cdf.(dvals); label = "inferred - $strat", line = linestyle(strat))
	end
	cdf = ecdf(pair_data[strategies[1]].distance_real)
	plot!(dvals, cdf.(dvals); label = "real", line = (:dash, :black))

	plot!(
		xlabel = "distance",
		ylabel = "",
		title = "Cumulative distribution of distance between leaves",
		legend = :bottomright,
		legendfontsize = 18,
		xlim = extrema(dvals) .* 1.1,
	)
end

# ╔═╡ 9167990a-33ed-41d9-a0c4-e3c6dfedf872
if fam_from_folder(folder_full) == fam_main
	fam = fam_from_folder(folder_full)
	method = method_from_folder(folder_full)

	p = plot(
		plt_dist_inf_v_real, plt_dist_cumulative;
		layout = grid(1,2),
		size = (1600, 800),
		margin = 10mm,
		dpi = 300,
	)

	if method == "arnet"
		savefig(p, joinpath(
			savedir, "branch_length_reconstruction_$(method)_$(fam).png"
		))
	end
end

# ╔═╡ Cell order:
# ╠═dea00e2e-8bea-4715-af4f-4b15ebeb792e
# ╠═81bc7c34-f1c5-11ee-17a1-f1983df9cbad
# ╠═ff4bc077-81bc-4dd7-80de-58ba09f4ff76
# ╠═4158b556-77a3-443c-96ef-2b90bab9ec97
# ╠═93396504-09f3-44a5-a23d-834784b41218
# ╠═4d06f815-17e9-41cf-9a5a-2c88e917665d
# ╠═94e3fc55-8a73-494d-8e5e-8f7aae8b554d
# ╠═8ed7f4ca-343f-45d5-8333-25ee836a81db
# ╠═7a763db5-9890-4c8d-b468-c93239a2aeb4
# ╠═c9981050-3a9c-4cda-b1f3-145de001713e
# ╟─6331fa1c-6d4c-43e2-94bd-d755f9a85e82
# ╠═53b74403-c8b2-4bb0-9038-50ef587bc2bf
# ╠═67bf7432-fd4d-410d-bb70-6e8344d7572e
# ╠═ecd42c56-6100-4a99-a044-011bc93245db
# ╠═6413baf8-deea-46e5-b609-896f95f7e4e8
# ╠═b93e7a60-9cbc-4f65-8083-ef58841ced47
# ╠═baee935b-9491-40d2-92eb-a648a190d1eb
# ╠═4a9f6f32-8ac5-43ee-8874-6a651546099f
# ╠═6ecb98ae-d792-4307-8a63-d212eb4c24fe
# ╟─b587e52b-5330-4ced-9a40-2cb39dac0ae4
# ╟─bd983086-a0ae-4072-b191-b1fcd63ed4e7
# ╠═1bc4405a-3249-4886-92e9-913b9516f221
# ╠═fc4dabcc-dfe7-411b-a485-2559ab136b5a
# ╠═26699dfd-04f5-40ce-a59d-142603b5c1fe
# ╟─9389795a-45fb-4ea5-b43d-c0599a6fe218
# ╠═ab30475e-a052-4dc4-b316-42b33d2d2205
# ╠═9167990a-33ed-41d9-a0c4-e3c6dfedf872
# ╟─6c407d72-2b03-4ae4-a95b-463fb438d711
# ╠═439679e1-b121-4f82-a68e-a311af5280ee
# ╠═2484add9-eb75-4a32-81e6-6a5ac68a6b39
# ╠═968c4424-fcc4-4d9e-b104-642e491e65f4
# ╠═1a2afd1a-6067-4a81-bda1-449c971d2060
# ╠═a86493bd-daac-4302-b378-9c7c65db7152
# ╠═825d2389-d53c-4d2c-8c87-e521f050c5e0
