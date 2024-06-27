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
	using JLD2
	using StatsPlots
	using StatsBase
end

# ╔═╡ 885ae0af-a181-40d8-b15c-ab8fc9f4d59d
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 186a0669-e014-4625-894f-2cf700006d03
let
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
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
	"Stiffler/subalignments", "seqstyle_uniform_treestyle_star_gencode"
)

# ╔═╡ e8bcceb8-5a57-45e0-bf7e-2dcf1d53b47e
data_nobin = let
	df = CSV.read(joinpath(folder_nobin, "data.csv"), DataFrame)
	# @info df
	gdf = @chain df @subset(@byrow !ismissing(:nerr_arnet)) groupby(:M)
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

# ╔═╡ e74d5941-8402-421e-941a-efc1088f55d3


# ╔═╡ f525b22e-194f-4dd6-ac6c-4c0732dc8d22


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
		# xlim = (1, 1000),
		ylim = (2, 16),
	)

	@df data_nobin scatter!(:M, :mean_arnet, label="arnet + genetic code")
	plot!(xlabel = "Number of leaves", ylabel = "Distance to w.t.")
end

# ╔═╡ afe621e4-9a32-4c7d-8b78-ca739d72b960
md"## Positions of errors"

# ╔═╡ 6700ff47-9ea1-49df-942d-ce1c84d8d829
begin
	strategies = [:cons, :iqtree, :arnet]
	strategy_names = ["consensus", "iqtree", "autoregressive"]
end

# ╔═╡ 98812576-010a-4648-aeb0-7f40d1305cc8


# ╔═╡ 03c4e063-83f0-4cf1-90aa-2a5182309932
md"## Likelihood"

# ╔═╡ 47aadd5c-636e-4c81-b822-0bb3a6477c58
begin
	lk_wt = readdlm(datadir("Stiffler/arnet/likelihood_wt.csv"))[1]
	lk_nat = readdlm(datadir("Stiffler/arnet/likelihood_nat.csv")) |> vec
	lk_R20 = readdlm(datadir("Stiffler/arnet/likelihood_R20.csv")) |> vec
end;

# ╔═╡ c8ce3ae0-84fe-42ca-8cf1-fd4f69d7c568
let p = plot()
	density!(lk_nat; line = (:black, 2), fill=(true, :black, .1), label = "")
	vline!([lk_wt], label="wild-type", color=1)
	density!(lk_R20, label="Round 20", color=2)
	plot!(
		xlim = (-350, -145),
		xlabel = "log-likelihood", 
	)
end

# ╔═╡ 123f98c2-c25a-4cef-9342-920954886f7d
md"# Utils"

# ╔═╡ 8186b193-114c-4a1a-9121-b514c1d92beb
begin
	string_to_vec(s::AbstractString) = readdlm(IOBuffer(s[2:end-1]), ',', Int) |> vec
	string_to_vec(x) = x
end

# ╔═╡ 2179ee13-c4e7-446d-b7d9-98e14cde160e
data_wcode = let
	# data contains results of arnet without the genetic code
	# this new df will replace this with arnet with the genetic code
	df_base = @chain "seqstyle_uniform_treestyle_star" begin
		datadir("Stiffler/subalignments", _)
		CSV.read(joinpath(_, "data.csv"), DataFrame)
		dropmissing
		@select(Not(r"arnet")) # remove arnet data
	end


	df_code = @chain "seqstyle_uniform_treestyle_star_gencode" begin
		datadir("Stiffler/subalignments", _)
		CSV.read(joinpath(_, "data.csv"), DataFrame)
		@subset(@byrow !ismissing(:nerr_arnet)) 
		@select(Not(r"iqtree"))
		@select(Not(r"cons"))
	end

	
	df = innerjoin(
		df_base, df_code, on = [:M, :id]; 
		makeunique=true, order = :left, renamecols = "" => ""
	)

	foreach(names(df)) do col_name
		if occursin(r"(rec)|(pos)", col_name)
			@transform!(df, $col_name = string_to_vec.($col_name))
		end
	end
	df
end

# ╔═╡ 6058a473-da5a-4a23-99d3-3dde0f951cc5
let
	readdir(datadir("Stiffler/subalignments"))
	outfolder = mkpath(datadir("Stiffler/subalignments/Results"))
	CSV.write(joinpath(outfolder, "data.csv"), data_wcode)
	@save(joinpath(outfolder, "data.jld2"), data_wcode)
	write(
		joinpath(outfolder, "info.txt"), 
		"""
			- tree is star
			- arnet is with gencode, and scaling branch length
			- iqtree uses JTT
		"""
	)
end

# ╔═╡ 7a0bba00-4c0b-4126-9d4b-a2bf7358dcd3
data_wcode_grouped = @chain data_wcode begin
	@subset(@byrow !ismissing(:nerr_arnet)) 
	groupby(:M)
	@combine(
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

# ╔═╡ b8a2cfe5-5529-4a2d-a048-57ca20c34c18
let p = plot()
	for strat in (:cons, :iqtree, :arnet)
		@df data_wcode_grouped plot!(:M, cols(Symbol(:lk_, strat)), label = "$strat")
	end
	hline!([lk_wt], label="wild-type", line = (:black, :dash, 3))
	plot!(
		xscale = :log10,
		xlabel = "M",
		ylabel = "log-likelihood",
	)

end

# ╔═╡ 5cf0f5a8-330f-4a6a-87fa-1aaa184a16ac
function get_error_positions(data, strat; M = 160)
	X = @chain data begin
		@subset :M .== M
		@select :X = cols(Symbol(:pos_, strat))
		@transform :X = map(string_to_vec, :X)
		flatten(:X)
		countmap(_.X)
		filter(x -> x[2] > 5, _)
	end
end

# ╔═╡ 15aba057-6997-4a0c-ad4d-dceaffec5538
error_positions, n_error_per_pos = let
	err_pos = Dict(
		strat => get_error_positions(data_wcode, strat) for strat in strategies
	)
	all_positions = mapreduce(vcat, strategies) do s 
		@chain err_pos[s] keys collect
	end |> unique |> sort
	error_per_pos = mapreduce(hcat, strategies) do s
		[get(err_pos[s], i, 0) for i in all_positions]
	end
	all_positions, error_per_pos
end

# ╔═╡ aac76244-7096-4aa3-b6c5-46f3cc511dfe
let p = plot()
	groupedbar(
		n_error_per_pos/100;
		bar_position = :dodge, 
		bar_width=0.5, 
		line = (2),
		xticks=(1:length(error_positions), error_positions),
		label = reshape(strategy_names, 1, 3),
	)

	plot!(
		ylim = (-0.01, 1.25),
		# xlim = (-.5, length(error_positions))
	)
end

# ╔═╡ Cell order:
# ╠═6a76c5ba-1ccd-11ef-1527-096e1c28dbf6
# ╠═885ae0af-a181-40d8-b15c-ab8fc9f4d59d
# ╠═186a0669-e014-4625-894f-2cf700006d03
# ╠═6f53d24c-d19a-48d5-8dbb-1bcde407d0a1
# ╠═fc88e1a6-5153-4428-9c1b-aac50f8ee9a0
# ╠═0fa2c8a3-e1e1-4d3d-8a6e-2d29905cd4ae
# ╠═040b87ca-9db5-4a44-9ce4-bda5426a1721
# ╠═e8bcceb8-5a57-45e0-bf7e-2dcf1d53b47e
# ╠═2179ee13-c4e7-446d-b7d9-98e14cde160e
# ╠═e74d5941-8402-421e-941a-efc1088f55d3
# ╠═6058a473-da5a-4a23-99d3-3dde0f951cc5
# ╠═f525b22e-194f-4dd6-ac6c-4c0732dc8d22
# ╠═7a0bba00-4c0b-4126-9d4b-a2bf7358dcd3
# ╟─50d03b54-0ceb-43de-881d-e737c6902bbe
# ╠═d7bd0229-1caf-4af2-9b8d-599f2ca8a1ea
# ╟─afe621e4-9a32-4c7d-8b78-ca739d72b960
# ╠═6700ff47-9ea1-49df-942d-ce1c84d8d829
# ╠═15aba057-6997-4a0c-ad4d-dceaffec5538
# ╠═5cf0f5a8-330f-4a6a-87fa-1aaa184a16ac
# ╠═aac76244-7096-4aa3-b6c5-46f3cc511dfe
# ╠═98812576-010a-4648-aeb0-7f40d1305cc8
# ╟─03c4e063-83f0-4cf1-90aa-2a5182309932
# ╟─47aadd5c-636e-4c81-b822-0bb3a6477c58
# ╠═b8a2cfe5-5529-4a2d-a048-57ca20c34c18
# ╠═c8ce3ae0-84fe-42ca-8cf1-fd4f69d7c568
# ╠═123f98c2-c25a-4cef-9342-920954886f7d
# ╠═8186b193-114c-4a1a-9121-b514c1d92beb
