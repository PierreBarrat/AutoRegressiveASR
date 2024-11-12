### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 4cac39d2-9dd0-11ef-14ba-f77e6111896a
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")
end

# ╔═╡ 86874393-2fc3-48b5-95aa-407f71b873b8
begin
	using AutoRegressiveASR
	using Chain
	using CSV
	using DataFrames
	using DataFramesMeta
	using Measures
	using PlutoUI
	using StatsBase
	using StatsPlots
end

# ╔═╡ 4ff8c740-d694-4cc1-a809-7d09fe104104
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 39799a1e-be2e-41f2-a237-e92a68e17d23
begin
	
	Plots.default(;pubfig(20)...)
end

# ╔═╡ d17b279d-5cec-44ea-b045-b6be3ce686ca
datdir = datadir("simulated/arnet_yule_rootposition")

# ╔═╡ 69df36ea-feb4-419e-8470-c0a87279ee60
folder_picker = @bind _folder Select(readdir(datdir))

# ╔═╡ 2dc36076-afd0-49a6-8fe5-84e66e72fbc4
folder = joinpath(datdir, _folder)

# ╔═╡ 277b92e0-7ac7-4247-b571-2e97a7bd7e07
data = let
	df = CSV.read(joinpath(folder, "measures.csv"), DataFrame)
	@transform! df begin
		:loss = :hamming_to_real - :original_hamming_to_real
	end
end

# ╔═╡ ff85d00a-7771-499b-96bc-6125011288ad
P2 = let
	# histogram of change in reconstruction for big root displacements (>1.5)
	# getting values
	vals = (@subset data :delta_root .> 1.5).hamming_to_original_reconstruction
	
	edges = collect(0.0:0.01:0.2)
	hvals = edges[2:end]
	H = fit(Histogram, vals, edges)
	# data for plotting - y are the normalized weights
	idx = findall(>(0), H.weights)
	x = hvals[idx]
	y = H.weights[idx] / sum(H.weights)

	for (i, frac) in enumerate(1 .- cumsum(y))
		frac = round(frac, sigdigits=2)
		@info "Fraction $frac have variation in Hamming of more than $(x[i])"
	end

	p = plot(x, y; label="", yscale=:log10)
	plot!(p;
		xlabel = "Hamming distance",
		ylabel = "",
		title = "Distribution of variation in reconstruction",
		xlim = (-0.0, 0.21),
		gridalpha = 0.05,
	)
end

# ╔═╡ 103a0879-72a7-492f-a697-466c3670e4e5
begin
	# smoothing alg
	w = 10
	outliers_right = 0.
	smoothing_alg = :hist
end

# ╔═╡ e4fa2000-8cea-4e1d-9d59-7239db04b06c
function sem(ystd, N)
	# error on mean calculated from sample standard deviation and number of samples
	# [-1.96, 196] has 95% of the mass for a normal distribution
	return 1.96 * ystd ./ sqrt.(N)
end

# ╔═╡ 0aec2ffa-b8d9-4842-96e7-8fd084619655
P1 = let p = plot()
	# Change in reconstruction
	x, y, ystd, N = ASRU.easy_smooth(
		data, :delta_root, :hamming_to_original_reconstruction; 
		w, alg=smoothing_alg, outliers_right, 
	)
	yerr = sem(ystd, N)
	plot!(x, y; ribbon = yerr, fillalpha=.2, label="Change in reconstruction",)
	
	# Loss in performance
	x, y, ystd, N = ASRU.easy_smooth(
		data, :delta_root, :loss; 
		w, alg=smoothing_alg, outliers_right, 
	)
	yerr = sem(ystd, N)
	plot!(x, y; ribbon = yerr, fillalpha=.2, label="Performance loss")

	plot!(
		xlim = (0,2),
		xlabel = "Distance between roots",
		ylabel = "Hamming distance",
		title = "Change in reconstruction vs. position of root",
		legend = :right,
		gridalpha = 0.05,
	)
end

# ╔═╡ a38ed47d-335c-4af5-a709-c634b77beb4a
panel = plot(
	P1, P2; 
	layout = grid(1,2),
	dpi = 300,
	size = (1800, 900),
	margin = 10mm,
)

# ╔═╡ 7f987e25-bef2-4e2c-975a-f9850059c479
let
	savename = projectdir("notes/article/figures/SI/root_position.pdf")
	savefig(panel, savename)
end

# ╔═╡ Cell order:
# ╠═4cac39d2-9dd0-11ef-14ba-f77e6111896a
# ╠═86874393-2fc3-48b5-95aa-407f71b873b8
# ╠═4ff8c740-d694-4cc1-a809-7d09fe104104
# ╠═39799a1e-be2e-41f2-a237-e92a68e17d23
# ╠═d17b279d-5cec-44ea-b045-b6be3ce686ca
# ╠═69df36ea-feb4-419e-8470-c0a87279ee60
# ╠═2dc36076-afd0-49a6-8fe5-84e66e72fbc4
# ╠═277b92e0-7ac7-4247-b571-2e97a7bd7e07
# ╟─0aec2ffa-b8d9-4842-96e7-8fd084619655
# ╟─ff85d00a-7771-499b-96bc-6125011288ad
# ╠═a38ed47d-335c-4af5-a709-c634b77beb4a
# ╠═7f987e25-bef2-4e2c-975a-f9850059c479
# ╠═103a0879-72a7-492f-a697-466c3670e4e5
# ╠═e4fa2000-8cea-4e1d-9d59-7239db04b06c
