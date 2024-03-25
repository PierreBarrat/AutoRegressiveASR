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

# ╔═╡ 988a4668-d718-11ee-1816-b73ece48c950
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using AutoRegressiveASR
	using DataFrames
	using DataFramesMeta
	using PlutoUI
	using StatsBase
	using StatsPlots
end

# ╔═╡ 3a06c22c-bbfe-41c1-9b68-26482fc1cc34
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 79db37a3-ed65-41f0-9e95-275bf7fac3a2
begin
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 13fef31f-abaa-4e94-917b-943b1b92901d
DIV = pluto_ingredients(scriptsdir("figures_and_results/diversity_functions.jl"))

# ╔═╡ 7990e53e-27cc-4d31-9d1b-c2c172d20d5c
_fs = @bind folder_full Select(readdir(datadir("simulated/potts_yule"); join=true))

# ╔═╡ 1e74de0b-656a-4341-a655-ebf8719de9fe
title_id = let
	folder = basename(folder_full)
	split(folder, "_")[1]
end

# ╔═╡ 7c3a2726-4925-4ed6-b8fc-dcc02f8f466d
path_folder, data_folder, prefix = let
	p, d = dirname(folder_full), basename(folder_full)
	prefix = split(d, "_")[1]
	p, d, prefix
end

# ╔═╡ 077bd775-60ce-41b3-bdaa-771b7784788d
data, filename = produce_or_load(
   Dict("basefolder" => folder_full, "out" => "diversity_data.jld2"); 
   filename = x -> joinpath(x["basefolder"], x["out"]),
   suffix = "",
) do config
   data = DIV.diversity_data(config["basefolder"], config["out"])
end;

# ╔═╡ 53a0574f-8def-49f3-b909-01465c5e8c26
md"# Figures"

# ╔═╡ 716af069-d904-4d5b-af49-d01c59949ec5
M = let
	M1 = @select data["iqtree"] @byrow :M = :Meff / :Meff_scaled
	M2 = @select data["ardca"] @byrow :M = :Meff / :Meff_scaled
	M = @chain begin 
		intersect(M1.M, M2.M)
		map(Int ∘ round, _)
		unique
	end
	length(M) > 1 && @warn "Is M well defined? $M"
	M[1]
end

# ╔═╡ de5239e6-c7b2-4e8f-a003-cd55fe3ef2c7
md"## Misc"

# ╔═╡ 4b1101bc-c54d-424e-848e-df40dad4feb8
strategies = ["iqtree", "ardca"]

# ╔═╡ 741a5584-8061-46ed-9922-ce0030651728
begin
	# smoothing width
	w = 20
	smoothing_alg = :hist
	pal = palette(:default)
	strat_color = Dict(strat => pal[i] for (i, strat) in enumerate(strategies))
end

# ╔═╡ 789e5665-173e-45af-8af2-2a194001fb95
names(data[strategies[1]])

# ╔═╡ 177dc4cb-b999-44fb-9600-e6ec3677c019
strat_label = Dict(
	"iqtree" => "iqtree",
	"ardca" => "autoregressive",
)

# ╔═╡ 600c4781-b372-4e69-a4d3-6ff27331bc58
linestyle = let
	lw = 4
	Dict(strat => (lw, strat_color[strat]) for strat in strategies)
end

# ╔═╡ e97ccac2-c5f1-4866-8d85-7f2a78dc691d
let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :depth, :Meff; w, alg=smoothing_alg,
		)
		plot!(x, y, label=strat_label[strat], line=linestyle[strat])
	end
	plot!(
		xlabel = "Depth",
		title = "$title_id, Meff (M=$M)",
	)
end

# ╔═╡ fe2bc5ee-1f28-4cdd-8521-1394f96f7faa
let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :depth, :av_self_hamming; w, alg=smoothing_alg,
		)
		plot!(x, y, label=strat_label[strat], line=linestyle[strat])
	end
	plot!(
		xlabel = "Depth",
		title = "$title_id - Average pw. hamming distance",
		# ylim = (0, .5),
	)
end

# ╔═╡ 6fcabb5d-d54b-4428-b04e-3158e6aecfb7
let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :depth, :entropy; w, alg=smoothing_alg,
		)
		plot!(x, y, label=strat_label[strat], line=linestyle[strat])
	end
	plot!(
		xlabel = "Depth",
		title = "$title_id - Entropy",
		# ylim = (0, .5),
	)
end

# ╔═╡ 58a81f2b-b728-4ca7-b8f8-9eeaa0cafb7f
markerstyle = let
	Dict(strat => (strat_color[strat], stroke(0)) for strat in strategies)
end

# ╔═╡ b27860a4-fde6-4903-befd-e039ade30d9d
let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :depth, :av_self_hamming; w, alg=smoothing_alg,
		)
		@df data[strat] scatter!(
			:depth, :av_self_hamming;
			label=strat_label[strat], marker=markerstyle[strat]
		)
		plot!(x, y, label="", line=linestyle[strat])
	end
	plot!(
		xlabel = "Depth",
		title = "$title_id - Average pw. hamming distance",
		# ylim = (-.05, .5),
		# xscale = :log10,
		legend = :topleft,
	)
end

# ╔═╡ 70a448e6-7e0f-4b4e-be30-05a5d4b30796
md"# Functions"

# ╔═╡ 8a4cab3c-892d-4d25-b72a-55f0f2da9be8
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

# ╔═╡ Cell order:
# ╠═988a4668-d718-11ee-1816-b73ece48c950
# ╠═3a06c22c-bbfe-41c1-9b68-26482fc1cc34
# ╠═79db37a3-ed65-41f0-9e95-275bf7fac3a2
# ╠═13fef31f-abaa-4e94-917b-943b1b92901d
# ╠═7990e53e-27cc-4d31-9d1b-c2c172d20d5c
# ╠═1e74de0b-656a-4341-a655-ebf8719de9fe
# ╠═7c3a2726-4925-4ed6-b8fc-dcc02f8f466d
# ╠═077bd775-60ce-41b3-bdaa-771b7784788d
# ╟─53a0574f-8def-49f3-b909-01465c5e8c26
# ╟─716af069-d904-4d5b-af49-d01c59949ec5
# ╟─e97ccac2-c5f1-4866-8d85-7f2a78dc691d
# ╠═741a5584-8061-46ed-9922-ce0030651728
# ╟─b27860a4-fde6-4903-befd-e039ade30d9d
# ╠═fe2bc5ee-1f28-4cdd-8521-1394f96f7faa
# ╠═6fcabb5d-d54b-4428-b04e-3158e6aecfb7
# ╠═789e5665-173e-45af-8af2-2a194001fb95
# ╟─de5239e6-c7b2-4e8f-a003-cd55fe3ef2c7
# ╠═4b1101bc-c54d-424e-848e-df40dad4feb8
# ╠═177dc4cb-b999-44fb-9600-e6ec3677c019
# ╠═600c4781-b372-4e69-a4d3-6ff27331bc58
# ╠═58a81f2b-b728-4ca7-b8f8-9eeaa0cafb7f
# ╠═70a448e6-7e0f-4b4e-be30-05a5d4b30796
# ╠═8a4cab3c-892d-4d25-b72a-55f0f2da9be8
