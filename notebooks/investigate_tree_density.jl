### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 2b69a640-fe5e-11ee-24e3-cd09d215792e
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using AutoRegressiveASR
	using BackwardCoalescent
	using CSV
	using DataFrames
	using DataFramesMeta
	using PlutoUI
	using StatsBase
	using StatsPlots
	using TreeTools
end

# ╔═╡ 77a6cc3c-072e-4edb-8448-063428dc2df5
heights = [.1, .25, .5, .75, 1., 1.5, 2]

# ╔═╡ 42ac4367-663c-497c-8343-2134004bd284
nleaves = 50

# ╔═╡ a05a2c57-a822-4377-8f91-679b2604d77e
function rescale_tree!(tree, H)
	T = TreeTools.distance(tree.root, first(leaves(tree)))
	foreach(n -> branch_length!(n, branch_length(n)*H/T), nodes(tree; skiproot=true))
	return nothing
end

# ╔═╡ 52af7224-dbb6-447b-874a-89f3e97f88e1
tree = genealogy(KingmanCoalescent(nleaves, 1))

# ╔═╡ a3bcb688-c15f-4555-a065-9b6a46cef41a
length(nodes(tree))

# ╔═╡ 406afede-1f2e-4f18-9ecc-f7dc2a4b3f68
τ = .1

# ╔═╡ d11442f3-ce14-4f29-8ab4-5935a9c01503
function depth_v_q(tree, τ)
	t = weighted_leaf_distance!(tree, τ)
	mapreduce(vcat, nodes(t)) do n 
		[n.data.depth n.data.Q]
	end
end

# ╔═╡ 3f2b6008-fe40-4e44-9392-2a5ae6091f08
kingman_dat = map(heights) do H
	mapreduce(vcat, 1:200) do rep
		tree = genealogy(KingmanCoalescent(nleaves, 1))
		rescale_tree!(tree, H)
		depth_v_q(tree, τ)
	end
end

# ╔═╡ 551f4e9b-b083-480d-8c0a-2373a97a3a2a
yule_dat = map(heights) do H
	mapreduce(vcat, 1:200) do rep
		tree = genealogy(YuleCoalescent(nleaves, 1))
		rescale_tree!(tree, H)
		depth_v_q(tree, τ)
	end
end

# ╔═╡ 5a22fb8b-73a8-42cb-8fe1-9715e421b41b
X = depth_v_q(tree, .1)

# ╔═╡ 19690666-2460-444b-9795-fac2e0d68f7d
distance(tree.root, first(leaves(tree)))

# ╔═╡ 57cf8706-f93d-472c-a107-3a7917c0f593
pal = palette(:bluesreds, length(heights))

# ╔═╡ 18e43af3-017f-403d-9655-48a551fd188b
begin
	# smoothing alg
	w = 20
	outliers_right = 0.
	smoothing_alg = :hist
end

# ╔═╡ 63715582-1791-4a76-8d55-d4c95b245df0
let p = plot()
	X = kingman_dat
	for (i, h) in enumerate(heights)
		x, y, _, _ = ASRU.easy_smooth(
			X[i][:, 1], X[i][:,2]; 
			w, alg=smoothing_alg, outliers_right
		)
		plot!(x, y, color = pal[i], label = "height = $h")
	end
	plot!(
		title = "Kingman", 
		xlabel = "depth",
		ylabel = "Information from leaves"
	)
end

# ╔═╡ 84bacdf4-810c-469d-a8c2-9e247154e555
let p = plot()
	X = yule_dat
	for (i, h) in enumerate(heights)
		x, y, _, _ = ASRU.easy_smooth(
			X[i][:, 1], X[i][:,2]; 
			w, alg=smoothing_alg, outliers_right
		)
		plot!(x, y, color = pal[i], label = "height = $h")
	end
	plot!(
		title = "Yule", 
		xlabel = "depth",
		ylabel = "Information from leaves"
	)
end

# ╔═╡ Cell order:
# ╠═2b69a640-fe5e-11ee-24e3-cd09d215792e
# ╠═77a6cc3c-072e-4edb-8448-063428dc2df5
# ╠═42ac4367-663c-497c-8343-2134004bd284
# ╠═a05a2c57-a822-4377-8f91-679b2604d77e
# ╠═52af7224-dbb6-447b-874a-89f3e97f88e1
# ╠═a3bcb688-c15f-4555-a065-9b6a46cef41a
# ╠═3f2b6008-fe40-4e44-9392-2a5ae6091f08
# ╠═551f4e9b-b083-480d-8c0a-2373a97a3a2a
# ╠═406afede-1f2e-4f18-9ecc-f7dc2a4b3f68
# ╠═d11442f3-ce14-4f29-8ab4-5935a9c01503
# ╠═5a22fb8b-73a8-42cb-8fe1-9715e421b41b
# ╠═19690666-2460-444b-9795-fac2e0d68f7d
# ╠═63715582-1791-4a76-8d55-d4c95b245df0
# ╠═84bacdf4-810c-469d-a8c2-9e247154e555
# ╠═57cf8706-f93d-472c-a107-3a7917c0f593
# ╠═18e43af3-017f-403d-9655-48a551fd188b
