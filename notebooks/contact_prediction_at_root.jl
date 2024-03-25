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

# ╔═╡ c1539d7a-e5de-11ee-3308-f97647f00678
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

# ╔═╡ 5f160aa8-2282-4606-8a83-10360803ddab
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 574892e8-92f9-4601-a6db-cd65a7783347
begin
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 7c650580-85d2-4307-84e5-b4f0a6592f41
_fs = @bind folder_full Select(readdir(datadir("simulated/potts_yule"); join=true))

# ╔═╡ 94e18f6d-3e2f-489e-b1ff-8384b30572ce


# ╔═╡ c7144fdf-57ec-4063-b9a5-1b6c8e32b514
CP = pluto_ingredients(scriptsdir("figures_and_results/contact_prediction.jl"))

# ╔═╡ 527c209d-9e0c-4127-8f2a-d831df601032
data, filename = produce_or_load(
	Dict("basefolder" => folder_full);
	filename = x -> joinpath(x["basefolder"], "contact_prediction.jld2"), 
	suffix = "",
) do config
	CP.contact_prediction_data(config["basefolder"])
end

# ╔═╡ 21347d90-cd04-4266-8dea-000f67cef5cc
data["ppv_eq"]

# ╔═╡ 7e4effea-8b76-4d18-a061-06929a082b1d
begin
	pal = palette(:default)
	color_iqtree = pal[1]
	color_ardca = pal[2]
	
	line_eq = (:dash, :black)
	line_iqtree = (color_iqtree)
	line_ardca = (color_ardca)
end;

# ╔═╡ d3c16917-3db9-4ab0-8d79-9593db37803e
let p = plot()
	plot!(data["ppv_eq"], line = line_eq, label = "")

	for (i, repdat) in data["asr"]
		# plot!(repdat["ardca"], line = line_iqtree, label = "")
	end
	Meff_ar = mean(x -> x["Meff_iqtree"], values(data["asr"])) |> round
	plot!(
		mean(x -> x["ppv_iqtree"], values(data["asr"]));
		line = line_iqtree, label = "iqtree - $(Meff_ar)",
	)

	Meff_ar = mean(x -> x["Meff_ardca"], values(data["asr"])) |> round
	plot!(
		mean(x -> x["ppv_ardca"], values(data["asr"]));
		line = line_ardca, label = "ArDCA - $(Meff_ar)",
	)
	
	plot!(xscale = :log10)
end

# ╔═╡ Cell order:
# ╠═c1539d7a-e5de-11ee-3308-f97647f00678
# ╠═5f160aa8-2282-4606-8a83-10360803ddab
# ╠═574892e8-92f9-4601-a6db-cd65a7783347
# ╠═7c650580-85d2-4307-84e5-b4f0a6592f41
# ╠═94e18f6d-3e2f-489e-b1ff-8384b30572ce
# ╠═c7144fdf-57ec-4063-b9a5-1b6c8e32b514
# ╠═527c209d-9e0c-4127-8f2a-d831df601032
# ╠═21347d90-cd04-4266-8dea-000f67cef5cc
# ╠═d3c16917-3db9-4ab0-8d79-9593db37803e
# ╠═7e4effea-8b76-4d18-a061-06929a082b1d
