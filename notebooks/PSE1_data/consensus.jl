### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 7f85291a-135e-11ef-23d2-95a813469efb
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__)

	using BioSequenceMappings
	using Chain
	using Distributions
	using Plots
	using StatsBase
end

# ╔═╡ f9e4d3d4-2df0-4a1c-96fc-a8f2bfc947e1
data_folder = datadir("Stiffler/aligned_data_ref")

# ╔═╡ a2702cca-49bc-409b-98fa-f426d9b6eb03
readdir(data_folder)

# ╔═╡ bc8a4f11-da89-4340-bb29-ae13f3e05cb3
alphabet = Alphabet(:aa);

# ╔═╡ 39f096f5-29ed-44c5-b627-253123525190
begin
	wt = read_fasta(joinpath(data_folder, "PSE1_aligned_PF13354_noinserts.fasta"))
	R20 = read_fasta(
		joinpath(data_folder, "PSE1_rnd20_aligned_PF13354_noinserts.fasta")
	)
	pfam = read_fasta(joinpath(data_folder, "PF13354_msa.fasta"))
end;

# ╔═╡ 1b4c5b12-868d-45b3-86bc-59db43e959e4
md"## Consensus of round 20"

# ╔═╡ ad804aa4-fe1a-417a-9151-9e33a9a1d923
f1_pfam = site_specific_frequencies(pfam; as_vec=false);

# ╔═╡ 893b19a4-9c8f-4b8a-83a3-42f8048c9a64
cons20 = consensus(R20)[1]

# ╔═╡ 01559685-f8f6-4f9f-ab3a-9735fa910ffe
positions = findall(cons20 .!= wt[1])

# ╔═╡ 08707927-d0f1-4da7-a7f1-f60cde1d7ad9
hamming(wt[1], cons20; normalize=false)

# ╔═╡ 68d1050f-05b6-44e3-925d-51d6aee9837c
map(positions) do pos
	aa_cons = cons20[pos]
	aa_wt = wt[1][pos]

	labels = ["." for a in 1:21]
	labels[aa_cons] = "cons"
	labels[aa_wt] = "wt"
	
	f = f1_pfam[:, pos]
	idx = sortperm(f; rev=true)
	f = f[idx]
	labels = labels[idx]
	hcat(labels, f)
end

# ╔═╡ 3d79ec84-ebc3-4c5e-ac96-c5f70a906f82
md"## Consensus of subalignments"

# ╔═╡ e291743e-53a6-4a53-9fb3-b407cb854673
Mvals = 1.5 .^ (5:15)

# ╔═╡ 7d3f8601-ef6e-4d83-a906-83487a7a4c2b
function hamming_Mcons(aln, M; nreps=100)
	return map(1:nreps) do _
		X = BioSequenceMappings.subsample_random(aln, M)
		hamming(wt[1], consensus(X)[1]; normalize=false)
	end
end

# ╔═╡ 7659ae77-1c4c-4c6e-a3e9-064a7901123b
dat_subaln = let
	Mvals = map(x -> Int(round(x)), 1.5 .^ (5:15))
	X = map(Mvals) do m  
		h = hamming_Mcons(R20, m)
		mean(h), std(h)
	end
	(Mvals, [x[1] for x in X], [x[2] for x in X])
end

# ╔═╡ 05ecd9ed-336c-477c-ab62-515b7215d805
let
	p = plot()
	scatter(dat_subaln[1], dat_subaln[2], label="")
	plot!(xscale = :log10)
end

# ╔═╡ bdec6623-2920-47e9-a69a-7633c223758e
md"## Hamming distance from round 20 to w.t."

# ╔═╡ e73568e4-dd6b-4bdf-9a70-ead02b70e8a7
L = size(wt, 1)

# ╔═╡ 8a8b2ca4-da99-4ba2-8d9a-26b5e0ccb085
H20 = pairwise_hamming(wt, R20; step_right=2) |> vec

# ╔═╡ 97e3a791-f324-449b-a78f-e6a1efe7363b
let p = plot()
	hvals = range(0, .25, length=500)
	# round 20
	X = ecdf(H20)
	plot!(hvals, X.(hvals), label = "Round 20")

	plot!(
		xlabel = "Hamming distance to w.t.",
	)

	# Poisson of mean <L*H20>
	P = Poisson(L*mean(H20))
	Hrand = ecdf(rand(P, 10000)/L)
	plot!(hvals, Hrand.(hvals), label = "poisson")
end

# ╔═╡ b3e5321d-f85d-4038-82a0-16d5163d76a1
let p = plot()
	hvals = 1:(L/4)
	# round 20
	X = ecdf(H20*L)
	plot!(hvals, X.(hvals), label = "Round 20")

	plot!(
		xlabel = "Hamming distance to w.t.",
	)
end

# ╔═╡ 23ddfc34-208f-4670-bb72-53f7f8638786
md"## Pairwise hamming distance of R20 sequences"

# ╔═╡ b891c8d0-7ddd-485c-af6f-a34f502dd2d4
md"""
"""

# ╔═╡ 70300524-18cf-4770-972f-73244ba6d2b1
τ = let
	H_wt_av = mean(H20) # average hamming to wild type
	-log(1 - H_wt_av) # expected branch length, crude approx.
end

# ╔═╡ c188e961-316d-4802-8fb3-583a28cae733
let p = plot()
	h = ecdf(pairwise_hamming(R20; step=100))
	hvals = range(0, .5, length=100)
	plot!(hvals, h.(hvals), label="data")

	# Poisson of parameter (1 - exp(-2τ))*L
	Hrand = @chain 1 - exp(-2*τ) begin
		_*L
		Poisson
		rand(10_000)
		_/L
		ecdf
	end
	plot!(hvals, Hrand.(hvals), label="weird theoretical curve")
end

# ╔═╡ Cell order:
# ╠═7f85291a-135e-11ef-23d2-95a813469efb
# ╠═f9e4d3d4-2df0-4a1c-96fc-a8f2bfc947e1
# ╠═a2702cca-49bc-409b-98fa-f426d9b6eb03
# ╠═bc8a4f11-da89-4340-bb29-ae13f3e05cb3
# ╠═39f096f5-29ed-44c5-b627-253123525190
# ╟─1b4c5b12-868d-45b3-86bc-59db43e959e4
# ╠═ad804aa4-fe1a-417a-9151-9e33a9a1d923
# ╠═893b19a4-9c8f-4b8a-83a3-42f8048c9a64
# ╠═01559685-f8f6-4f9f-ab3a-9735fa910ffe
# ╠═08707927-d0f1-4da7-a7f1-f60cde1d7ad9
# ╠═68d1050f-05b6-44e3-925d-51d6aee9837c
# ╟─3d79ec84-ebc3-4c5e-ac96-c5f70a906f82
# ╠═7659ae77-1c4c-4c6e-a3e9-064a7901123b
# ╠═e291743e-53a6-4a53-9fb3-b407cb854673
# ╠═05ecd9ed-336c-477c-ab62-515b7215d805
# ╠═7d3f8601-ef6e-4d83-a906-83487a7a4c2b
# ╟─bdec6623-2920-47e9-a69a-7633c223758e
# ╠═e73568e4-dd6b-4bdf-9a70-ead02b70e8a7
# ╠═8a8b2ca4-da99-4ba2-8d9a-26b5e0ccb085
# ╠═97e3a791-f324-449b-a78f-e6a1efe7363b
# ╟─b3e5321d-f85d-4038-82a0-16d5163d76a1
# ╟─23ddfc34-208f-4670-bb72-53f7f8638786
# ╠═b891c8d0-7ddd-485c-af6f-a34f502dd2d4
# ╠═70300524-18cf-4770-972f-73244ba6d2b1
# ╠═c188e961-316d-4802-8fb3-583a28cae733
