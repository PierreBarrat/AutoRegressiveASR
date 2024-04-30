### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ cd7fdf66-06cd-11ef-1102-ebc33fe74b51
begin
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using BackwardCoalescent
	using BenchmarkTools
	using Plots
	using StatsBase
	using TreeTools
end

# ╔═╡ eec251ae-235a-4462-8f94-b17f4f1135ca
begin
	local plot_defaults = joinpath(homedir(), ".julia/config/plot_defaults.jl")
	if isfile(plot_defaults)
		include(plot_defaults)
	end
end

# ╔═╡ d2948310-b60f-4074-8666-b55f3262ed4b
md"# Setup"

# ╔═╡ a39c54bd-6687-4438-a8d6-ec40f2178743
if isdefined(@__MODULE__, :pubfig)
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 28a44a02-082f-41e8-989d-6f2c7af5a03c
begin
	figdir = mkpath(joinpath(plotsdir(), "paper/SI/"))
	figdir_paper = joinpath(projectdir(), "notes/article/figures/SI")
end

# ╔═╡ df085193-9fbd-434a-bede-2eeb477529be
md"# Figures"

# ╔═╡ 9e029eae-f178-4adb-bf74-787051d81dc7
begin
	n = 100 # number of leaves in trees
	Nreps = 1000 # number of repetitions of each tree
end

# ╔═╡ 776bf117-d633-4e42-84b0-ba7778bb423c
md"# Functions"

# ╔═╡ e4d5d759-606f-4342-8d4c-4024f8452551
function fit_histogram(edges, D::AbstractVector)
    normalize(x) = x / sum(x)
    W = fit(Histogram, D, edges).weights |> normalize
    xvals = (edges[2:end] + edges[1:end-1]) / 2
    return xvals, W
end

# ╔═╡ 6249659f-7c43-41ef-80a7-7129cb6c5b6c
function _append_depths!(D, node::TreeNode, h, i)
	if isleaf(node)
		return i
	end
	
	h += ismissing(branch_length(node)) ? 0 : branch_length(node)
	i += 1
	D[i] = h
	for c in children(node)
		i = _append_depths!(D, c, h, i)
	end
	return i
end

# ╔═╡ 71c7ac1b-a96c-4801-9704-d92e7343e543
function rescale_tree!(tree)
	H = distance(root(tree), first(leaves(tree)))
	foreach(nodes(tree)) do node
		τ = branch_length(node)
		!ismissing(τ) && branch_length!(node, τ/H)
	end
end

# ╔═╡ ae8c0f1d-6a02-4a1d-a306-090f5a671332
function get_node_depths!(tree)
	rescale_tree!(tree)
	H = distance(root(tree), first(leaves(tree)))
	N = length(nodes(tree)) - length(leaves(tree))
	D = Vector{Float64}(undef, N)
	_append_depths!(D, root(tree), 0, 0)
	return H .- D
end

# ╔═╡ 68596bff-0893-41a8-87dd-f66c52803554
node_depths_kingman = let
	coa = KingmanCoalescent(n, 1)
	mapreduce(vcat, 1:Nreps) do _ 
		tree = genealogy(coa)
		get_node_depths!(tree)
	end
end	

# ╔═╡ abaee924-00e5-4ea7-873e-abb51a2aa065
node_depths_yule = let
	coa = YuleCoalescent(n, 1)
	mapreduce(vcat, 1:Nreps) do _ 
		tree = genealogy(coa)
		get_node_depths!(tree)
	end
end	

# ╔═╡ 1714bd40-e978-4763-a9b0-cb1dd9118e52
plt = let p = plot()
	edges = 0:.025:1
	dvals, depth_distribution_kingman = fit_histogram(edges, node_depths_kingman)
	plot!(
		dvals, depth_distribution_kingman;
		label = "Kingman", 
	)

	dvals, depth_distribution_yule = fit_histogram(edges, node_depths_yule)
	plot!(
		dvals, depth_distribution_yule;
		label = "Yule", 
	)

	plot!(
		yscale = :log10, 
		ylim = (1e-4, 1),
		xlabel = "Node depth",
	)

	savefig(joinpath(figdir_paper, "depth_distribution_coalescents.png"))
	p
end

# ╔═╡ 98fb5846-0975-49b0-82db-7afcdcf200bc
let p = plot()
	kg = ecdf(node_depths_kingman)
	yule = ecdf(node_depths_yule)

	dvals = 0:.01:1
	plot!(dvals, kg.(dvals))
	plot!(dvals, yule.(dvals))
end

# ╔═╡ f9039b52-9e3f-40e9-97e3-72e8b20261c8
let
	coalescent = KingmanCoalescent(100, 1)
	tree = genealogy(coalescent)
	X1 = get_node_depths!(tree)
	X2 = map(internals(tree)) do node
		TreeTools.distance_to_closest_leaf(tree, label(node))
	end
	scatter(sort(X1), sort(X2))
	plot!([0.,1], [0,1], line=(:black, :dash))
	title!("Check that my function works")
end

# ╔═╡ 518a7a28-4963-4f78-8d69-622600fa582c
let
	coalescent = KingmanCoalescent(100, 1)
	tree = genealogy(coalescent)
	@time X1 = get_node_depths!(tree)
end

# ╔═╡ Cell order:
# ╠═cd7fdf66-06cd-11ef-1102-ebc33fe74b51
# ╠═d2948310-b60f-4074-8666-b55f3262ed4b
# ╠═eec251ae-235a-4462-8f94-b17f4f1135ca
# ╠═a39c54bd-6687-4438-a8d6-ec40f2178743
# ╠═28a44a02-082f-41e8-989d-6f2c7af5a03c
# ╟─df085193-9fbd-434a-bede-2eeb477529be
# ╠═9e029eae-f178-4adb-bf74-787051d81dc7
# ╠═68596bff-0893-41a8-87dd-f66c52803554
# ╠═abaee924-00e5-4ea7-873e-abb51a2aa065
# ╠═1714bd40-e978-4763-a9b0-cb1dd9118e52
# ╠═98fb5846-0975-49b0-82db-7afcdcf200bc
# ╟─776bf117-d633-4e42-84b0-ba7778bb423c
# ╠═e4d5d759-606f-4342-8d4c-4024f8452551
# ╠═f9039b52-9e3f-40e9-97e3-72e8b20261c8
# ╠═518a7a28-4963-4f78-8d69-622600fa582c
# ╠═ae8c0f1d-6a02-4a1d-a306-090f5a671332
# ╠═6249659f-7c43-41ef-80a7-7129cb6c5b6c
# ╠═71c7ac1b-a96c-4801-9704-d92e7343e543
