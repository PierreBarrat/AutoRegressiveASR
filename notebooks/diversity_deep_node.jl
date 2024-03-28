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

# ╔═╡ 3964dddc-eb71-11ee-25c8-7b3f4ffb63dd
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using AutoRegressiveASR
	using AncestralSequenceReconstruction
	using Chain
	using DataFrames
	using DataFramesMeta
	using DCATools
	using PlutoUI
	using StatsBase
	using StatsPlots
end

# ╔═╡ 61668a7d-b174-40bc-81bb-7c7c0b30ae36
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 33fcb755-77a6-4127-9cda-856d23db58d8
begin
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ ea32ae08-4dc8-49aa-b7df-21639d8b1732
target_node = "internal_48"

# ╔═╡ ead00618-12d9-4c91-9426-80822e7cd2d3


# ╔═╡ be6f4068-6795-4c7d-857b-87cc8453b13a
_fs = @bind folder_base Select(readdir(datadir("simulated/potts_yule"); join=true))

# ╔═╡ 9123275e-769f-4483-a961-eceaf28e1661
md"## Precomputed"

# ╔═╡ 4057f7d2-bb7f-4a4e-9573-1e88091c7e26
md"# Figures"

# ╔═╡ 992f9d52-d79a-41f8-a897-4a22adec3925
folder = let
	rep = 7
	joinpath(folder_base, "data/$rep")
end

# ╔═╡ 321a278f-6969-4e0f-ad35-78b80eb85211
alignment_leaves = read_msa(joinpath(folder, "alignment_leaves.fasta"));

# ╔═╡ d7af4857-d52e-4dd5-99b5-0715da2c8f27
seq_real = let
	aln = read_msa(joinpath(folder, "alignment_internals.fasta"))
	aln[target_node] |> collect
end

# ╔═╡ 50a92ab1-d7b9-4535-a184-02223952af84
iqtree_profile = @chain begin
	joinpath(folder, "iqtree/IQTREE.state")
	ASRU.parse_iqtree_state_file
	getindex(_, target_node)
	getproperty(_, :model)
end

# ╔═╡ d37e8e2f-ca06-4efe-a617-68ec463230fb


# ╔═╡ 4d87eac0-2ae7-4124-aa2d-5dfcc6b953f1
md"# Misc"

# ╔═╡ 858affc5-23f0-468a-abb8-3fdabfc9616f
begin
	pal = palette(:default)
	line_iqtree = (pal[1], 3)
	line_ar = (pal[2], 3)
	line_leaves = (:black, :dash)
end

# ╔═╡ a9519f81-2099-4b3e-bb0c-401f3ae2e230


# ╔═╡ 4302a755-dfe7-462a-8c38-850b38075e76
ASRU.hamming

# ╔═╡ fb00aa23-78f3-42c2-98d9-befc9bfdaaeb
function alignment_from_profile(profile::ASR.ProfileModel; M = 500)
    return @chain begin
        ASR.sample(profile, M)
        map(x -> ASR.intvec_to_sequence(x, profile.alphabet), eachcol(_))
        map(x -> ASR.sequence_to_intvec(x, ASR.Alphabet(:aa)), _)
        hcat(_...)'
        DCASample(_; mapping = ASR.reverse_mapping(ASR.Alphabet(:aa)))
    end
end

# ╔═╡ 06889df5-91fa-4cf0-b984-27ce47d8af4c
# Bayes sequences
begin
	alignment_ar = read_msa(joinpath(
		folder, "autoregressive_diversity", "Bayes", "reconstructed_$(target_node).fasta"
	))
	alignment_iqtree = alignment_from_profile(iqtree_profile)
end;

# ╔═╡ 682b265d-724a-492e-b903-1fe24e1ead47
begin
	pw_hamming_real = DCATools.pw_hamming_distance(
		alignment_leaves; exclude_gaps=true
	) |> ecdf
	
	pw_hamming_iqtree = DCATools.pw_hamming_distance(
		alignment_iqtree; exclude_gaps=true
	) |> ecdf

	pw_hamming_ar = DCATools.pw_hamming_distance(
		alignment_ar; exclude_gaps=true
	) |> ecdf
end

# ╔═╡ 9ed060de-e094-45fb-b998-77f1652f6499
plt_hamming_to_real = let p = plot()
	iqtree = map(alignment_iqtree) do x 
		DCATools.hamming(seq_real, x; exclude_gaps=true, normalize=true)
	end |> ecdf
	ar = map(alignment_ar) do x 
		DCATools.hamming(seq_real, x; exclude_gaps=true, normalize=true)
	end |> ecdf
	leaves = map(alignment_leaves) do x 
		DCATools.hamming(seq_real, x; exclude_gaps=true, normalize=true)
	end |> ecdf

	xvals = range(0, 1, 100)
	plot!(xvals, iqtree.(xvals), line = line_iqtree, label="iqtree")
	plot!(xvals, ar.(xvals), line = line_ar, label="autoregressive")
	plot!(xvals, leaves.(xvals), line = line_leaves, label="")

	plot!(
		xlabel = "Hamming distance to real",
		ylabel = "",
		title = "",
	)
end

# ╔═╡ 58d1e299-2395-43c2-af44-47673b57ca67
DCATools.pw_hamming_distance(alignment_iqtree; exclude_gaps=true) |> ecdf

# ╔═╡ 04bd2a62-5f3c-47d8-adec-5e822efa3aac
function ml_from_profile(profile::ASR.ProfileModel)
    return @chain begin ASR.ml_sequence(profile)
        ASR.intvec_to_sequence(profile.alphabet)
        ASR.sequence_to_intvec(ASR.Alphabet(:aa))
    end
end

# ╔═╡ 5b20a485-e654-404e-b516-deaf5bd99557
# ML sequences
begin
	mlseq_iqtree = ml_from_profile(iqtree_profile)
	mlseq_ar = @chain target_node begin
		joinpath(folder, "autoregressive_diversity/ML/reconstructed_$(_).fasta")
		read_msa
		first
		collect
	end
end;

# ╔═╡ a89e6fb9-c89e-4536-a01f-75a5d2dbeb0a
let p = plot()
	iqtree = map(alignment_iqtree) do x 
		DCATools.hamming(mlseq_iqtree, x; exclude_gaps=true, normalize=true)
	end |> ecdf
	ar = map(alignment_ar) do x 
		DCATools.hamming(mlseq_ar, x; exclude_gaps=true, normalize=true)
	end |> ecdf

	xvals = range(0, 1, 100)
	plot!(xvals, iqtree.(xvals), line = line_iqtree, label="iqtree")
	plot!(xvals, ar.(xvals), line = line_ar, label="autoregressive")

	plot!(
		xlabel = "Hamming distance to ML",
		ylabel = "",
		title = "",
	)
end

# ╔═╡ Cell order:
# ╠═3964dddc-eb71-11ee-25c8-7b3f4ffb63dd
# ╠═61668a7d-b174-40bc-81bb-7c7c0b30ae36
# ╠═33fcb755-77a6-4127-9cda-856d23db58d8
# ╠═ea32ae08-4dc8-49aa-b7df-21639d8b1732
# ╠═ead00618-12d9-4c91-9426-80822e7cd2d3
# ╠═be6f4068-6795-4c7d-857b-87cc8453b13a
# ╠═321a278f-6969-4e0f-ad35-78b80eb85211
# ╠═d7af4857-d52e-4dd5-99b5-0715da2c8f27
# ╠═50a92ab1-d7b9-4535-a184-02223952af84
# ╠═06889df5-91fa-4cf0-b984-27ce47d8af4c
# ╠═5b20a485-e654-404e-b516-deaf5bd99557
# ╟─9123275e-769f-4483-a961-eceaf28e1661
# ╠═682b265d-724a-492e-b903-1fe24e1ead47
# ╟─4057f7d2-bb7f-4a4e-9573-1e88091c7e26
# ╠═992f9d52-d79a-41f8-a897-4a22adec3925
# ╟─9ed060de-e094-45fb-b998-77f1652f6499
# ╠═a89e6fb9-c89e-4536-a01f-75a5d2dbeb0a
# ╠═d37e8e2f-ca06-4efe-a617-68ec463230fb
# ╠═4d87eac0-2ae7-4124-aa2d-5dfcc6b953f1
# ╠═858affc5-23f0-468a-abb8-3fdabfc9616f
# ╠═58d1e299-2395-43c2-af44-47673b57ca67
# ╠═a9519f81-2099-4b3e-bb0c-401f3ae2e230
# ╠═4302a755-dfe7-462a-8c38-850b38075e76
# ╠═fb00aa23-78f3-42c2-98d9-befc9bfdaaeb
# ╠═04bd2a62-5f3c-47d8-adec-5e822efa3aac
