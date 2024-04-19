### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 6f04b439-eb29-4969-819b-6e37233a13ec
begin
	using Revise
	using DrWatson
	@quickactivate("AutoRegressiveASR")

	using AutoRegressiveASR
	using BioSequenceMappings
	using Chain
	using DataFrames
	using DelimitedFiles
	using JLD2
	using Plots
	using Tulip
	using StatsBase
	using StatsPlots
end

# ╔═╡ d75eee42-5eb3-41bb-bd85-b440d2152b1c
basefolder = datadir("dynamics/PF00014_M=250_Ninit=100")

# ╔═╡ 6512a2b1-5b01-4fbf-91e2-5cc21f9a4a47
function eval_chain(
	basefolder, tvals, s0, aln_ref, self_hamming_ref, cross_hamming_ref,
)
	D = mapreduce(vcat, tvals) do t
		# X ~ sequences at time t
		X = read_fasta(joinpath(basefolder, "alignment_t$(t).fasta"))
		hcat(
			AutoRegressiveASR.emd_self_hamming(X, self_hamming_ref), 
			AutoRegressiveASR.emd_cross_hamming(X, aln_ref, cross_hamming_ref),
			mean(s -> hamming(s, s0), X),
		)
	end
	return D
end

# ╔═╡ 8287bdac-c8b3-463c-bb56-6cbd09aa1a1e
function split_ref_alignment(X; M_max = 500)
	M = length(X)
	m = min(M_max, div(M, 2))

	X1 = subsample_random(X, m)
	X2 = subsample_random(X, m)
	return X1, X2
end

# ╔═╡ 688c8f4c-7773-49d5-92b9-8a5b12e87348
function dynamics_measures(basefolder)
	params = JLD2.load(joinpath(basefolder, "parameters.jld2"))
	
	aln_nat = read_fasta(joinpath(basefolder, "alignment_nat.fasta"))
	aln_ref_1, aln_ref_2 = split_ref_alignment(aln_nat)
	self_hamming_ref = pairwise_hamming(aln_ref_1; as_vec=true)
	cross_hamming_ref = pairwise_hamming(aln_ref_1, aln_ref_2) |> vec

	initial_sequences = read_fasta(joinpath(basefolder, "init_sequences.fasta"))
	Nreps = length(initial_sequences)
	
	data = Dict()
	models = ["arnet", "potts"]
	for model in models 
		chain = model * "_chains"
		tvals = @chain joinpath(basefolder, chain, "") begin
			readdir(; join=true)
			filter(f -> occursin("time_values", f), _)
			first
			readdlm
			@. occursin("potts", chain) ? Int(_) : _
		end
		D = mean(1:2) do i
			eval_chain(
				joinpath(basefolder, chain, "$i"),
				tvals,
				initial_sequences[i],
				aln_ref_1,
				self_hamming_ref,
				cross_hamming_ref,
			)
		end
		data[model] = hcat(tvals, D)
	end

	dataframes = map(models) do model
		d = data[model]
		DataFrame(
			:t => d[:, 1],
			:emd_self_hamming => d[:, 2],
			:emd_cross_hamming => d[:, 3],
			:hamming_to_init => d[:, 4],
		)
	end
	return dataframes
end

# ╔═╡ 01ae6bd1-18b1-4408-af0d-d5311d294df8
data = dynamics_measures(basefolder)

# ╔═╡ 58adbb87-cd3f-470a-aa82-2e970f62996a
let p = plot()
	# @df data[1] plot!(:t, :emd_self_hamming)
	@df data[2] plot!(:t, :hamming_to_init)
	# plot!(yscale=:log)#, xlim = (1e-3, 2000))
end

# ╔═╡ 037bfc64-a5d8-4d1e-aa68-143ce6bcb58c
# ╠═╡ disabled = true
#=╠═╡
data = ans
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═6f04b439-eb29-4969-819b-6e37233a13ec
# ╠═d75eee42-5eb3-41bb-bd85-b440d2152b1c
# ╠═688c8f4c-7773-49d5-92b9-8a5b12e87348
# ╠═6512a2b1-5b01-4fbf-91e2-5cc21f9a4a47
# ╠═8287bdac-c8b3-463c-bb56-6cbd09aa1a1e
# ╠═01ae6bd1-18b1-4408-af0d-d5311d294df8
# ╠═58adbb87-cd3f-470a-aa82-2e970f62996a
# ╠═037bfc64-a5d8-4d1e-aa68-143ce6bcb58c
