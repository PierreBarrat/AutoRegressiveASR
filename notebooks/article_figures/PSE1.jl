### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 669af142-2ee7-11ef-3a61-059bd4388f68
begin
	using DrWatson
	quickactivate(@__DIR__)

	using Chain
	using CSV
	using DataFrames
	using DataFramesMeta
	using JLD2
	using Measures
	using StatsBase
	using StatsPlots
end

# ╔═╡ 6c2d97a1-144d-4744-9477-5dc524d66dde
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 2bab1daa-faeb-46aa-9b0a-67e4ad375242
let
	plt_defaults = pubfig(22)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 4022bae7-dfe3-4685-a210-3f85e35bb1ec
data = let
	@load datadir("Stiffler/subalignments/Results/data.jld2") data_wcode
	dropmissing(data_wcode)
end;

# ╔═╡ 5910cb98-fc83-41aa-9c55-80732b395e19
methods = (:cons, :iqtree, :arnet)

# ╔═╡ f369f07e-2db1-4c7d-962c-e59705b4d3fc
begin
	method_name = Dict(
		:cons => "consensus",
		:iqtree => "IQ-TREE",
		:arnet => "autoregressive",
	)

	# colors
	pal = palette(:default)
	method_color = Dict(
		:cons => pal[3],
		:iqtree => pal[1],
		:arnet => pal[2],
	)
end

# ╔═╡ d4c9337c-bca9-44ed-93bf-1a8595053436
md"# Panel"

# ╔═╡ ebb25fc3-86e3-41f6-9951-f762fe154b1c


# ╔═╡ d49c24dc-1c23-40a1-9810-10f1506ffec1
md"# Number of errors vs M"

# ╔═╡ d59c4ccf-b189-4de4-82cb-30fc119c9277
data_grouped = @chain data begin
	groupby(:M)
	@combine(
		:mean_cons = mean(:nerr_cons),
		:std_cons = std(:nerr_cons),
		:lk_cons = mean(:likelihood_cons),
		
		:mean_arnet = mean(:nerr_arnet),
		:std_arnet = std(:nerr_arnet),
		:lk_arnet = mean(:likelihood_arnet),
		
		:mean_iqtree = mean(:nerr_iqtree),
		:std_iqtree = std(:nerr_iqtree),
		:lk_iqtree = mean(:likelihood_iqtree),
		
		:N = length(:nerr_cons),
	)
end;

# ╔═╡ f67a63e3-6e1b-4755-8ea7-c56c5a76aa8a
plt_H_v_M = let p = plot()
	for m in methods
		@df data_grouped scatter!(
			:M, cols(Symbol(:mean_, m));
			marker = (:o, 10, stroke(.5, :black)),
			color =  method_color[m],
			yerr = (cols(Symbol(:std_, m)) ./ sqrt.(:N)),
			label = method_name[m],
		)
	end
	plot!(
		xscale = :log10,
		ylim = (0, 15),
		xlabel = "Number of leaves",
		ylabel = "Hamming distance",
		xticks = ([10, 100], ["10", "100"])
	)
end

# ╔═╡ 3ecb35ad-bd8e-40ce-9454-d068ab14caec
data_grouped

# ╔═╡ 716347c5-93c3-43f2-b9b5-6696c17ea267
md"# Errors per position"

# ╔═╡ 91d05a7f-91b6-46d7-a04b-b76d0d0f3b55
Mref = 640

# ╔═╡ ca81e407-2db5-49a1-9343-a6d633fe5176
function get_error_positions(data, strat, Mref)
	X = @chain data begin
		@subset :M .== Mref
		@select :X = cols(Symbol(:pos_, strat))
		# @transform :X = map(string_to_vec, :X)
		flatten(:X)
		countmap(_.X)
		filter(x -> x[2] > 10, _)
	end
end

# ╔═╡ 50c3dab6-5186-4605-abd8-7101bf7ddecc
error_positions, n_error_per_pos = let
	err_pos = Dict(
		m => get_error_positions(data, m, Mref) for m in methods
	)
	all_positions = mapreduce(vcat, methods) do s 
		@chain err_pos[s] keys collect
	end |> unique |> sort
	error_per_pos = mapreduce(hcat, methods) do s
		[get(err_pos[s], i, 0) for i in all_positions]
	end
	all_positions, error_per_pos
end

# ╔═╡ 7c9329b4-18b4-4c2d-bd95-e2fb1ad47336
plt_errors_per_position = let p = plot()
	groupedbar(
		n_error_per_pos/100;
		bar_position = :dodge, 
		bar_width=0.5, 
		line = (2),
		xticks=(1:length(error_positions), error_positions),
		label = reshape([method_name[m] for m in methods], 1, 3),
		color = reshape([method_color[m] for m in methods], 1, 3)
	)

	plot!(
		ylim = (-0.01, 1.3),
		yticks = range(0, 1., 6),
		xlabel = "Position",
		ylabel = "Fraction of errors",
	)
end

# ╔═╡ a3e060eb-9449-4408-8608-f87688ffed2d
panel = plot(
	plt_H_v_M, plt_errors_per_position,
	layout = grid(1, 2),
	size = (1600, 800), 
	bottom_margin = 15mm, left_margin = 10mm,
	dpi = 300,
)


# ╔═╡ 72f71c1f-84ab-4ffd-93e6-509d5c8d9692
begin
	savename = "stiffler_reconstruction_left.png"
	savedir = projectdir("notes/article/figures/")
	savefig(panel, projectdir(savedir, savename))
end

# ╔═╡ Cell order:
# ╠═669af142-2ee7-11ef-3a61-059bd4388f68
# ╠═6c2d97a1-144d-4744-9477-5dc524d66dde
# ╠═2bab1daa-faeb-46aa-9b0a-67e4ad375242
# ╠═4022bae7-dfe3-4685-a210-3f85e35bb1ec
# ╠═5910cb98-fc83-41aa-9c55-80732b395e19
# ╠═f369f07e-2db1-4c7d-962c-e59705b4d3fc
# ╠═d4c9337c-bca9-44ed-93bf-1a8595053436
# ╠═a3e060eb-9449-4408-8608-f87688ffed2d
# ╠═ebb25fc3-86e3-41f6-9951-f762fe154b1c
# ╠═72f71c1f-84ab-4ffd-93e6-509d5c8d9692
# ╟─d49c24dc-1c23-40a1-9810-10f1506ffec1
# ╠═d59c4ccf-b189-4de4-82cb-30fc119c9277
# ╠═f67a63e3-6e1b-4755-8ea7-c56c5a76aa8a
# ╠═3ecb35ad-bd8e-40ce-9454-d068ab14caec
# ╟─716347c5-93c3-43f2-b9b5-6696c17ea267
# ╠═91d05a7f-91b6-46d7-a04b-b76d0d0f3b55
# ╠═ca81e407-2db5-49a1-9343-a6d633fe5176
# ╠═50c3dab6-5186-4605-abd8-7101bf7ddecc
# ╠═7c9329b4-18b4-4c2d-bd95-e2fb1ad47336
