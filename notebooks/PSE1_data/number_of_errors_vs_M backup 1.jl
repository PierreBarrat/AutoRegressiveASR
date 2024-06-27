### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 6a76c5ba-1ccd-11ef-1527-096e1c28dbf6
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__)

	using Chain
	using CSV
	using DataFrames
	using DataFramesMeta
	using DelimitedFiles
	using Distributions
	using StatsPlots
	using StatsBase
end

# ╔═╡ 6f53d24c-d19a-48d5-8dbb-1bcde407d0a1
base_folder = datadir("Stiffler/subalignments/seqstyle_uniform_treestyle_star")

# ╔═╡ fc88e1a6-5153-4428-9c1b-aac50f8ee9a0
df = CSV.read(joinpath(base_folder, "data.csv"), DataFrame) |> dropmissing;

# ╔═╡ 0fa2c8a3-e1e1-4d3d-8a6e-2d29905cd4ae
data = let
	gdf = groupby(df, :M)
	@combine(
		gdf, 
		:mean_cons = mean(:nerr_cons),
		:std_cons = std(:nerr_cons),
		:lk_cons = mean(:likelihood_cons),
		
		:mean_arnet = mean(:nerr_arnet),
		:std_arnet = std(:nerr_arnet),
		:lk_arnet = mean(:likelihood_arnet),
		
		:mean_iqtree = mean(:nerr_iqtree),
		:std_iqtree = std(:nerr_arnet),
		:lk_iqtree = mean(:likelihood_iqtree),
		
		:N = length(:nerr_cons),
	)
end;

# ╔═╡ 040b87ca-9db5-4a44-9ce4-bda5426a1721
folder_nobin = datadir(
	"Stiffler/subalignments", "seqstyle_uniform_treestyle_star_modelfinder"
)

# ╔═╡ e8bcceb8-5a57-45e0-bf7e-2dcf1d53b47e
data_nobin = let
	df = CSV.read(joinpath(folder_nobin, "data.csv"), DataFrame)
	# @info df
	gdf = @chain df @subset(@byrow !ismissing(:nerr_iqtree)) groupby(:M)
	@combine(
		gdf, 
		:mean_cons = mean(:nerr_cons),
		:std_cons = std(:nerr_cons),
		:lk_cons = mean(:likelihood_cons),
		
		:mean_arnet = mean(:nerr_arnet),
		:std_arnet = std(:nerr_arnet),
		:lk_arnet = mean(:likelihood_arnet),
		
		:mean_iqtree = mean(:nerr_iqtree),
		:std_iqtree = std(:nerr_arnet),
		:lk_iqtree = mean(:likelihood_iqtree),
		
		:N = length(:nerr_cons),
	)
end

# ╔═╡ 50d03b54-0ceb-43de-881d-e737c6902bbe
md"## Hamming to wt"

# ╔═╡ d7bd0229-1caf-4af2-9b8d-599f2ca8a1ea
let p = plot()
	for label in (:cons, :iqtree, :arnet)
		@df data scatter!(
			:M, cols(Symbol(:mean_, label));
			yerr = (cols(Symbol(:std_, label)) ./ sqrt.(:N)), label = "$label",
		)
	end
	plot!(
		xscale = :log10,
		xlim = (1, 1000),
		ylim = (3, 16),
	)

	@df data_nobin scatter!(:M, :mean_arnet)
end

# ╔═╡ 03c4e063-83f0-4cf1-90aa-2a5182309932
md"## Likelihood"

# ╔═╡ 47aadd5c-636e-4c81-b822-0bb3a6477c58
begin
	lk_wt = readdlm(datadir("Stiffler/arnet/likelihood_wt.csv"))[1]
	lk_nat = readdlm(datadir("Stiffler/arnet/likelihood_nat.csv")) |> vec
	lk_R20 = readdlm(datadir("Stiffler/arnet/likelihood_R20.csv")) |> vec
end;

# ╔═╡ b8a2cfe5-5529-4a2d-a048-57ca20c34c18
let p = plot()
	for strat in (:cons, :iqtree, :arnet)
		@df data plot!(:M, cols(Symbol(:lk_, strat)), label = "$strat")
	end
	hline!([lk_wt], label="wild-type", line = (:black, :dash))
	plot!(
		xscale = :log10,
		xlabel = "M",
		ylabel = "log-likelihood",
	)

end

# ╔═╡ b9ea548d-2550-44ef-8b80-2143748eea65
# side plot with lk distribution
density_lk_eq = let
	lk_lim = (-270, -200)
	plt = density(lk_nat)#; xlim=lk_lim)
	# x, y = plt[1][1][:x], plt[1][1][:y]
	# plot(
	# 	y, x;
	# 	ylim=lk_lim, label="", xticks = [], yticks = [], line=(3, :black),
	# 	frame=:no, xaxis=false, fill=true, fillcolor = :black, fillalpha=.1
	# )
end

# ╔═╡ 12108eb5-e469-4724-a3a2-88880715721b
StatsBase.std([1.,2.])

# ╔═╡ Cell order:
# ╠═6a76c5ba-1ccd-11ef-1527-096e1c28dbf6
# ╠═6f53d24c-d19a-48d5-8dbb-1bcde407d0a1
# ╠═fc88e1a6-5153-4428-9c1b-aac50f8ee9a0
# ╠═0fa2c8a3-e1e1-4d3d-8a6e-2d29905cd4ae
# ╠═040b87ca-9db5-4a44-9ce4-bda5426a1721
# ╠═e8bcceb8-5a57-45e0-bf7e-2dcf1d53b47e
# ╟─50d03b54-0ceb-43de-881d-e737c6902bbe
# ╠═d7bd0229-1caf-4af2-9b8d-599f2ca8a1ea
# ╟─03c4e063-83f0-4cf1-90aa-2a5182309932
# ╠═47aadd5c-636e-4c81-b822-0bb3a6477c58
# ╠═b8a2cfe5-5529-4a2d-a048-57ca20c34c18
# ╠═b9ea548d-2550-44ef-8b80-2143748eea65
# ╠═12108eb5-e469-4724-a3a2-88880715721b
