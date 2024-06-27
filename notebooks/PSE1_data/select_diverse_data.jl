### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ c1a44edc-1c07-11ef-2a01-9d06ecdbad4a
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__)

	using BioSequenceMappings
	using Chain
	using Distributions
	using Plots
	using Random
	using StatsBase
end

# ╔═╡ 596b7e62-7080-44c5-88f8-531ed2b25e6e
data_folder = datadir("Stiffler/aligned_data_ref")

# ╔═╡ 5ce30c41-d11c-4f43-865d-d36d166cb2a1
	R20 = read_fasta(
		joinpath(data_folder, "PSE1_rnd20_aligned_PF13354_noinserts.fasta")
	);

# ╔═╡ 3e69c9c4-948b-4fc9-a12f-11c732d0e739


# ╔═╡ 2ee5e5a9-fe84-41d7-9f35-0d995d48dfef
function sample_diverse!(X, aln, threshold)
	idx = randperm(length(aln))
	for s in eachcol(aln[idx])
		min_dist = minimum(x -> hamming(x, s), X)
		if min_dist > threshold
			push!(X, collect(s))
			return true
		end
	end
	return false
end

# ╔═╡ dfcf6fe3-97c4-43fc-ae71-dfcf256c9f0d
function sample_diverse(aln, M; threshold = 0.1)
	X = map(collect, rand(aln, 1))
	for m in 1:M-1 
		flag = sample_diverse!(X, aln, threshold)
		if !flag
			@warn "Stopped at $(m) sequences"
			break
		end
	end
	return Alignment(X, Alphabet(aln))
end

# ╔═╡ a3916d85-7f20-439a-b467-06d1b4bf580d
X = sample_diverse(R20, 100; threshold = .2)

# ╔═╡ bee65442-ccd6-4ef7-8de7-a7a6e82de468
pairwise_hamming(X)

# ╔═╡ Cell order:
# ╠═c1a44edc-1c07-11ef-2a01-9d06ecdbad4a
# ╠═596b7e62-7080-44c5-88f8-531ed2b25e6e
# ╠═5ce30c41-d11c-4f43-865d-d36d166cb2a1
# ╠═3e69c9c4-948b-4fc9-a12f-11c732d0e739
# ╠═dfcf6fe3-97c4-43fc-ae71-dfcf256c9f0d
# ╠═2ee5e5a9-fe84-41d7-9f35-0d995d48dfef
# ╠═a3916d85-7f20-439a-b467-06d1b4bf580d
# ╠═bee65442-ccd6-4ef7-8de7-a7a6e82de468
