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

# ╔═╡ a8ec28f8-c9ad-11ee-1efc-751443eb65d6
begin
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")
	using AutoRegressiveASR
	using Chain
	using DCATools
	using JSON3
	using Plots
	using PlutoUI
	using StatsBase
	using TreeTools
end

# ╔═╡ c1391eae-6dd7-4b34-a3d9-fdaa3a1f4830
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 2b763d1c-0863-4792-be6e-fc705e1be322
let
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ adc5c0a6-c599-4f85-bfa3-6654b4578478
md"# Pairwise hamming distance for simulated data"

# ╔═╡ 026ba3fd-65dd-4d4d-bdc1-79761d17b47c
folder_list = vcat(
	readdir(datadir("simulated/arnet_yule"); join=true),
	readdir(datadir("simulated/potts_yule"); join=true),
)

# ╔═╡ 3b4abe01-a079-4ca8-95e7-94be09ecac7a
folder_picker = @bind _tmp_dat Select(folder_list)

# ╔═╡ 79495f33-e493-4400-9435-725793a058d3
dat_folder = datadir("simulated/arnet_yule", _tmp_dat)

# ╔═╡ da0768e4-41e7-4404-8a36-14d5978ca4dd
potts, sample_eq = let
	parameters = JSON3.read(
		open(joinpath(dat_folder, "simulation_parameters.json"), "r"), 
	)
	(
		nothing, #DCAGraph(parameters[:generative_model]), 
		read_msa(parameters[:sample_equilibrium])
	)
end

# ╔═╡ e755a5d7-ac63-4fe1-a9a6-b399073294f3
pw_hamming_eq = DCATools.pw_hamming_distance(sample_eq; step=5) |> ecdf

# ╔═╡ 40d86f59-414f-4b2a-8898-48eb6e87f660
figdir = mkpath(joinpath(plotsdir(), "potts_yule", basename(dat_folder)))

# ╔═╡ 257ea632-a1ad-40e1-818a-a04906bc1170
pw_hamming_trees = map(ASRU.get_tree_folders(joinpath(dat_folder, "data"))) do fol
	aln = joinpath(fol, "alignment_leaves.fasta") |> read_msa
	DCATools.pw_hamming_distance(aln) |> ecdf
end

# ╔═╡ b618a701-c11b-456d-b747-c7161afb879c
let p = plot()
	xvals = 0:.01:1
	plot!(xvals, pw_hamming_eq.(xvals), label="Eq.", line=(:black))
	for d in pw_hamming_trees
		plot!(xvals, d.(xvals), color=1, label="", alpha = .5, lw = 1)
	end
	plot!(
		xlabel = "Hamming distance",
		title = "Distribution of pairwise hamming distances",
	)

	savefig(p, joinpath(figdir, "pw_hamming_aln_leaves.png"))
	p
end

# ╔═╡ a43ef7a1-cc46-4e0c-8e5d-93455d78a1f7
md"## Closest ancestral sequence for each leaf"

# ╔═╡ df59e5b6-8968-4758-b39c-5d3c6f0e1a2f
closest_eq_sequence = map(DCATools.subsample(sample_eq, 1:10:10_000)) do s
	minimum(DCATools.subsample(sample_eq, 2:10:10_000)) do x 
		DCATools.hamming(s, x)/length(x)
	end
end |> ecdf

# ╔═╡ 2ddb0d89-81d4-46af-a102-c4a2826823ac
closest_ancestral_seq = let
	X = map(ASRU.get_tree_folders(joinpath(dat_folder, "data"))) do fol
		leaves = joinpath(fol, "alignment_leaves.fasta") |> read_msa
		internals = joinpath(fol, "alignment_internals.fasta") |> read_msa
		map(leaves) do leaf
			minimum(x -> DCATools.hamming(x, leaf), internals)/length(leaf)
		end
	end
	vcat(X...)
end |> ecdf

# ╔═╡ 9d129f47-b77d-4395-9e07-d9dba60416ba
let p = plot()
	xvals = 0:.01:1
	plot!(xvals, closest_eq_sequence.(xvals), label="Equilibrated sample", line=(:black))
	plot!(xvals, closest_ancestral_seq.(xvals), label="")
	plot!(
		xlabel = "Hamming",
		title = "CDF: Hamming distance to closest internal",
	)
end

# ╔═╡ 4c8eadbb-2a4b-4da4-bd80-696b58d4015e
	parameters = JSON3.read(
		open(joinpath(dat_folder, "simulation_parameters.json"), "r"), 
	)

# ╔═╡ a930c0ac-03db-4cc9-a670-a3881601f41b
# ╠═╡ disabled = true
#=╠═╡
	parameters = JSON3.read(
		open(joinpath(dat_folder, "simulation_parameters.json"), "r"), 
	) |> keys |> collect
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═a8ec28f8-c9ad-11ee-1efc-751443eb65d6
# ╠═c1391eae-6dd7-4b34-a3d9-fdaa3a1f4830
# ╠═2b763d1c-0863-4792-be6e-fc705e1be322
# ╠═79495f33-e493-4400-9435-725793a058d3
# ╠═da0768e4-41e7-4404-8a36-14d5978ca4dd
# ╠═4c8eadbb-2a4b-4da4-bd80-696b58d4015e
# ╠═40d86f59-414f-4b2a-8898-48eb6e87f660
# ╠═a930c0ac-03db-4cc9-a670-a3881601f41b
# ╟─adc5c0a6-c599-4f85-bfa3-6654b4578478
# ╠═e755a5d7-ac63-4fe1-a9a6-b399073294f3
# ╠═257ea632-a1ad-40e1-818a-a04906bc1170
# ╠═026ba3fd-65dd-4d4d-bdc1-79761d17b47c
# ╠═3b4abe01-a079-4ca8-95e7-94be09ecac7a
# ╠═b618a701-c11b-456d-b747-c7161afb879c
# ╟─a43ef7a1-cc46-4e0c-8e5d-93455d78a1f7
# ╠═df59e5b6-8968-4758-b39c-5d3c6f0e1a2f
# ╠═2ddb0d89-81d4-46af-a102-c4a2826823ac
# ╠═9d129f47-b77d-4395-9e07-d9dba60416ba
