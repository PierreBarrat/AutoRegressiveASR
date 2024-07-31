### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 8001ce72-4f39-11ef-0e06-2d89322f54f0
begin
	using DrWatson
	quickactivate(@__DIR__)

	using ArDCA
	using BioSequenceMappings
	using Chain
	using CSV
	using DataFrames
	using DataFramesMeta
	using JLD2
	using StatsBase
	using StatsPlots
end

# ╔═╡ e75ab66b-3a32-488a-ba21-e6f17e4bf90e
arnet = JLD2.load(datadir(
	"Stiffler/arnet/arnet_PF13354_lJ0.01_lH0.001.jld2"
))["arnet"]

# ╔═╡ b9873f2a-c99c-4e6c-a2fc-221a09505020
data = let
	@load datadir("Stiffler/subalignments/Results/data.jld2") data_wcode
	dropmissing(data_wcode)
end;

# ╔═╡ de143533-ca60-4311-a5e1-778a619470aa
wt_to_swissprot = let
	dat = CSV.read(
		datadir("Stiffler/subalignments/Results/list_muts.remap.tsv"), 
		DataFrame
	)
	Dict(r.matteo_pierre => r.PSE1_swissprot for r in eachrow(dat))
end

# ╔═╡ 2a08a98b-d502-415a-90e9-1cb95aa7ee2f
methods = (:cons, :iqtree, :arnet)

# ╔═╡ c995c8d4-1602-49f9-9cc5-9178c5c63266
Mref = 640

# ╔═╡ 7920bd1c-8946-4675-8430-5a797ad800cd
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

# ╔═╡ 9c3fa3af-36c4-4f03-a0c3-5e46ee292030
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

# ╔═╡ 7e58a830-c1ed-4286-914b-c8b6a9dc736c
wt_to_uniprot(i) = wt_to_swissprot[i]

# ╔═╡ f2577725-a4a2-419e-be1f-9444ca2025de
md"## Table of mutations in wt context"

# ╔═╡ 17aa7b6f-2b66-47e2-8866-3a1f21e3468a
mutations = map(error_positions) do i
	mutations = []
	for m in methods
		subdat = @chain data begin
			@subset :M .== Mref
			@select begin
				:err = cols(Symbol(:rec_aa_, m))
				:pos = cols(Symbol(:pos_, m))
			end
		end
		for r in eachrow(subdat)
			idx = findfirst(==(i), r.pos)
			!isnothing(idx) && push!(mutations, r.err[idx])
		end
	end
	countmap(mutations)
	@chain mutations countmap findall(>(10), _)
end

# ╔═╡ 37269b0c-af45-4f0e-8582-d84aafb297a7
wt = read_fasta(datadir(
	"Stiffler/aligned_data_ref/PSE1_aligned_PF13354_noinserts.fasta"
))[1] |> copy

# ╔═╡ 5ef98dc6-dddd-4eca-b66d-490cf8513db6
wt_state = [wt[i] for i in error_positions]

# ╔═╡ 75be7478-a196-43c7-b735-03110a7205c0
alphabet = Alphabet(:aa)

# ╔═╡ ac1ad2b3-33d9-40de-8f4f-570350ba8d9a
dat_mut = let
	df = DataFrame(pos=Int[], wt=Char[], mut=Char[], delta_loglk = Float64[])
	loglk_wt = ArDCA.loglikelihood(wt, arnet)
	for (i, m) in zip(error_positions, mutations)
		mut = copy(wt)
		mut[i] = m[1]
		r = Dict(
			:pos => i,
			:wt => alphabet(wt[i]),
			:mut => alphabet(m[1]),
			:delta_loglk => ArDCA.loglikelihood(mut, arnet) - loglk_wt
		)
		push!(df, r)
	end
	df
end

# ╔═╡ 0744365a-2621-4bc8-92f1-5af3f7578ba4
pwd()

# ╔═╡ 46fae5dd-6dc2-41e9-969c-4c0d43b8f0de
md"## DMS and observed mutations (Δlk)"

# ╔═╡ 0b81c3d4-b351-4e54-ba0d-90296e9e26e6
function dms(wt, arnet)
	L = length(wt)
	q = 21
	lk_wt = ArDCA.loglikelihood(wt, arnet)
	
	X = zeros(Union{Missing,Float64}, q, L)
	m = copy(wt)
	for i in 1:L, a in 1:q
		m[i] = a
		X[a, i] = a == wt[i] ? missing : ArDCA.loglikelihood(m, arnet) - lk_wt
		m[i] = wt[i]
	end
	return X
end

# ╔═╡ 099a021a-b8f4-488e-8178-bc693b62ee98
dms(wt, arnet)[:, 212]

# ╔═╡ a8ed08b2-2a3b-46f4-b809-a819da1b52f2
alphabet(18)

# ╔═╡ bc96a636-59be-4441-bcfb-aae36cbe3a4f
let p = plot()
	dms_vals = @chain reshape(dms(wt, arnet), 1, :) vec skipmissing collect
	density!(dms_vals, label="")
	for r in eachrow(dat_mut)
		vline!([r.delta_loglk], label=r.pos)
	end
	p
end

# ╔═╡ 6768c7bd-b56d-436c-9cea-d430313b0640
md"## C --> S mutation: bad in w.t., reverse in a R20 seq.? "

# ╔═╡ ef34f74c-b69e-42bc-b42c-69e8f3c78ed1
rnd20 = read_fasta(datadir(
	"Stiffler/aligned_data_ref/PSE1_rnd20_aligned_PF13354_noinserts.fasta"
))

# ╔═╡ 79d2e22f-54b5-4528-b575-9844b6c8f082
site_specific_frequencies(rnd20; as_vec=false)[:, 212]

# ╔═╡ d60f41a8-6444-41bc-8ecc-7c0b2a985e2a
dms_C_to_S = let
	# Conclusion: putting back the serine is good! 
	i = 212
	mutants = findall(x -> alphabet(x[212]) == 'C', rnd20)
	map(mutants[1:50:end]) do m
		rnd20_mutant = rnd20[m]
		seq = copy(rnd20_mutant)
		seq[i] = alphabet('S')
		ArDCA.loglikelihood(seq, arnet) - ArDCA.loglikelihood(rnd20_mutant, arnet)
	end
end

# ╔═╡ Cell order:
# ╠═8001ce72-4f39-11ef-0e06-2d89322f54f0
# ╠═e75ab66b-3a32-488a-ba21-e6f17e4bf90e
# ╠═b9873f2a-c99c-4e6c-a2fc-221a09505020
# ╠═de143533-ca60-4311-a5e1-778a619470aa
# ╠═2a08a98b-d502-415a-90e9-1cb95aa7ee2f
# ╠═c995c8d4-1602-49f9-9cc5-9178c5c63266
# ╠═7920bd1c-8946-4675-8430-5a797ad800cd
# ╠═9c3fa3af-36c4-4f03-a0c3-5e46ee292030
# ╠═7e58a830-c1ed-4286-914b-c8b6a9dc736c
# ╟─f2577725-a4a2-419e-be1f-9444ca2025de
# ╠═17aa7b6f-2b66-47e2-8866-3a1f21e3468a
# ╠═37269b0c-af45-4f0e-8582-d84aafb297a7
# ╠═5ef98dc6-dddd-4eca-b66d-490cf8513db6
# ╠═75be7478-a196-43c7-b735-03110a7205c0
# ╠═ac1ad2b3-33d9-40de-8f4f-570350ba8d9a
# ╠═0744365a-2621-4bc8-92f1-5af3f7578ba4
# ╟─46fae5dd-6dc2-41e9-969c-4c0d43b8f0de
# ╠═0b81c3d4-b351-4e54-ba0d-90296e9e26e6
# ╠═099a021a-b8f4-488e-8178-bc693b62ee98
# ╠═a8ed08b2-2a3b-46f4-b809-a819da1b52f2
# ╠═bc96a636-59be-4441-bcfb-aae36cbe3a4f
# ╠═6768c7bd-b56d-436c-9cea-d430313b0640
# ╠═ef34f74c-b69e-42bc-b42c-69e8f3c78ed1
# ╠═79d2e22f-54b5-4528-b575-9844b6c8f082
# ╠═d60f41a8-6444-41bc-8ecc-7c0b2a985e2a
