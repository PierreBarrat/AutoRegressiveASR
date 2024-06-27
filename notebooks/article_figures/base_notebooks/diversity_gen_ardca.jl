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

# ╔═╡ 6c6ca152-0858-11ef-2581-a998a9aef0a8
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

# ╔═╡ 4e79b620-78b9-4a3e-83fc-1e84a887399d
begin
	local plot_defaults = joinpath(homedir(), ".julia/config/plot_defaults.jl")
	if isfile(plot_defaults)
		include(plot_defaults)
	end
end

# ╔═╡ e39b3831-f5f7-4a45-a92c-f704aa07c207
PlutoUI.TableOfContents()

# ╔═╡ 39b29c1e-1009-47f6-8a3c-1c37465f17bd
if isdefined(@__MODULE__, :pubfig)
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 05bda6a8-b577-4dd3-b193-db895b99a863
LM = pluto_ingredients(
	scriptsdir("figures_and_results/analyze_results_and_write_df.jl")
)

# ╔═╡ 7339904f-8b69-4dd2-9b4d-0bb5ffeba9ff
md"# Reading data"

# ╔═╡ 39bd18d7-1751-4acc-9852-f4eedf0c220a
md"## Folders"

# ╔═╡ d6359377-93dc-48ff-b14d-cb06e2d55b7b
folder_list = vcat(
	readdir(datadir("simulated/arnet_yule"); join=true),
)

# ╔═╡ d9ac189a-9ecd-40ee-abd5-8bb75042cab6
folder_picker = @bind folder_full Select(folder_list)

# ╔═╡ 183a61a3-879f-4431-a23e-6f18df4ee5cf
folder = basename(folder_full)

# ╔═╡ 02f26dce-1679-4d7d-b11b-bd589492c750
begin
	figdir = mkpath(joinpath(
		plotsdir(), "paper/diversity/ardca_yule/", basename(folder_full)
	))
	mkpath(figdir)
end

# ╔═╡ a03fd7ba-3cb5-4fbe-a419-a28179653958
md"## Reading data"

# ╔═╡ 0ef95b17-50fc-49e4-94c8-6d793eceb69b
DIV = pluto_ingredients(scriptsdir("figures_and_results/diversity_functions.jl"))

# ╔═╡ 9d58522c-a74f-4ffd-879c-780ca551b618
title_id = let
	folder = basename(folder_full)
	split(folder, "_")[1]
end

# ╔═╡ c6416a96-6e9b-4a1a-b3d3-7b2e4b2be37e
path_folder, data_folder, prefix = let
	p, d = dirname(folder_full), basename(folder_full)
	prefix = split(d, "_")[1]
	p, d, prefix
end

# ╔═╡ 64b1d353-8ec2-4338-bb5a-beb1d781f87a
data, filename = produce_or_load(
   Dict("basefolder" => folder_full, "out" => "diversity_data.jld2"); 
   filename = x -> joinpath(x["basefolder"], x["out"]),
   suffix = "",
) do config
   data = DIV.diversity_data(config["basefolder"], config["out"])
end;

# ╔═╡ 6c568f66-22df-4d7c-8810-1456babc6d11
md"# Figures"

# ╔═╡ 33ad28c7-0d7e-484d-9dd1-3edef25ed11f
function sem(ystd, N)
    # error on mean calculated from sample standard deviation and number of samples
    # [-1.96, 196] has 95% of the mass for a normal distribution
    return 1.96 * ystd ./ sqrt.(N)
end

# ╔═╡ dbc9b16d-0b7b-4908-801e-3fa2783dfec0
md"## Misc"

# ╔═╡ 12c33299-9c3d-4d9e-8969-42d707fb641f
strategies = ["iqtree", "ardca"]

# ╔═╡ fe521dba-36c1-4aad-a2f2-771c0a1a7aed
begin
	# smoothing width
	w = 20
	outliers_right = 0.
	smoothing_alg = :hist
	pal = palette(:default)
	strat_color = Dict(strat => pal[i] for (i, strat) in enumerate(strategies))
end

# ╔═╡ e99f52db-67dc-4a58-98d5-4d8f2563ca57
strat_label = Dict(
	"iqtree" => "iqtree",
	"ardca" => "autoregressive",
)

# ╔═╡ 98781a05-cab8-4f9f-9905-948c4c00c806
linestyle = let
	lw = 4
	Dict(strat => (lw, strat_color[strat]) for strat in strategies)
end

# ╔═╡ f0ca8ece-850a-427c-9b7a-07133a3ed96b
let p = plot()
	for strat in strategies
		x, y, ystd, N = ASRU.easy_smooth(
			data[strat], :depth, :av_self_hamming; 
			w, alg=smoothing_alg, outliers_right,
		)
		yerr = sem(ystd, N)
		plot!(
			x, y, ribbon = yerr;
			fillalpha = .2, label=strat_label[strat], line=linestyle[strat]
		)
	end
	plot!(
		xlabel = "Node depth",
		xlim = (-0.025, 2.025),
		ylabel = "Self-Hamming distance",
		title = "",
		frame = :box,
		legend = :topleft,
	)
	savefig(joinpath(figdir, "pwdistance_v_depth.png"))
	p
end

# ╔═╡ 8cd823ed-016a-4832-a2ea-83dc496b5603
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
	)
	savefig(joinpath(figdir, "entropy_v_depth.png"))
	p
end

# ╔═╡ 58f63e1e-41cf-444e-bfeb-ef2f1dc8596f
markerstyle = let
	Dict(strat => (strat_color[strat], stroke(0)) for strat in strategies)
end

# ╔═╡ Cell order:
# ╠═e39b3831-f5f7-4a45-a92c-f704aa07c207
# ╠═6c6ca152-0858-11ef-2581-a998a9aef0a8
# ╠═4e79b620-78b9-4a3e-83fc-1e84a887399d
# ╠═39b29c1e-1009-47f6-8a3c-1c37465f17bd
# ╠═05bda6a8-b577-4dd3-b193-db895b99a863
# ╠═7339904f-8b69-4dd2-9b4d-0bb5ffeba9ff
# ╠═39bd18d7-1751-4acc-9852-f4eedf0c220a
# ╠═d6359377-93dc-48ff-b14d-cb06e2d55b7b
# ╠═183a61a3-879f-4431-a23e-6f18df4ee5cf
# ╠═d9ac189a-9ecd-40ee-abd5-8bb75042cab6
# ╠═02f26dce-1679-4d7d-b11b-bd589492c750
# ╟─a03fd7ba-3cb5-4fbe-a419-a28179653958
# ╠═0ef95b17-50fc-49e4-94c8-6d793eceb69b
# ╠═9d58522c-a74f-4ffd-879c-780ca551b618
# ╠═c6416a96-6e9b-4a1a-b3d3-7b2e4b2be37e
# ╠═64b1d353-8ec2-4338-bb5a-beb1d781f87a
# ╟─6c568f66-22df-4d7c-8810-1456babc6d11
# ╠═33ad28c7-0d7e-484d-9dd1-3edef25ed11f
# ╠═f0ca8ece-850a-427c-9b7a-07133a3ed96b
# ╠═8cd823ed-016a-4832-a2ea-83dc496b5603
# ╠═dbc9b16d-0b7b-4908-801e-3fa2783dfec0
# ╠═fe521dba-36c1-4aad-a2f2-771c0a1a7aed
# ╠═12c33299-9c3d-4d9e-8969-42d707fb641f
# ╠═e99f52db-67dc-4a58-98d5-4d8f2563ca57
# ╠═98781a05-cab8-4f9f-9905-948c4c00c806
# ╠═58f63e1e-41cf-444e-bfeb-ef2f1dc8596f
