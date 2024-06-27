### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ d920b090-24a5-11ef-0ef8-8bd20e780e52
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__)

	using ArDCA
	using BioSequenceMappings
	using Chain
	using CSV
	using DataFrames
	using DataFramesMeta
	using DelimitedFiles
	using JLD2
	using StatsPlots
	using StatsBase
end

# ╔═╡ 4547b9d5-74b6-4d9e-8cab-6a9919ae3155
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 2100eba7-5edb-41b2-be3c-9465afb3cec2
let
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 815f6c9c-921c-40fa-bc90-4a8b26e79cdb
base_folder = datadir("Stiffler/subalignments/Results")

# ╔═╡ fdce6423-651f-4e3b-8737-f4191439c5b1
begin
	M_ref = 160
	alphabet = Alphabet(:aa)
end

# ╔═╡ 624397f5-9159-4482-a401-8cf71d1ddc55
df = @chain CSV.read(joinpath(base_folder, "data.csv"), DataFrame) dropmissing;

# ╔═╡ 6c76669d-2c05-4c5e-bb8d-648278e91c4c
begin
	strategies = [:cons, :iqtree, :arnet]
	strategy_names = ["consensus", "iqtree", "autoregressive"]
end

# ╔═╡ eaf5c62b-df07-45e1-8147-9581f36a66d5
arnet = JLD2.load(
	datadir("Stiffler/arnet/arnet_PF13354_lJ0.01_lH0.001.jld2")
)["arnet"]

# ╔═╡ 34ee72a4-6e8f-4479-835b-c1cd6a5d6dfb
md"## Work"

# ╔═╡ fa9e3a3e-afc9-44f6-b09e-739442243d44
md"# Functions & tools"

# ╔═╡ 9ff151f8-6807-4fef-8127-073f34a41285
function ar_local_field(s, i, arnet)
	j = findfirst(==(i), arnet.idxperm) - 1
	if j == 0
		return arnet.p0
	end
	
	H = copy(arnet.H[j])
	for k in 1:(j-1)
		H += arnet.J[j][:, s[arnet.idxperm[k]], k]
	end
	return ArDCA.softmax(H)
end

# ╔═╡ ad72228a-3381-4eef-9c15-75afa42694d0
normalize(X) = X / sum(X)

# ╔═╡ eef3db46-fce9-45bd-a906-67fa233cdce6
begin
	local folder = datadir("Stiffler/aligned_data_ref")
	wt = read_fasta(joinpath(folder, "PSE1_aligned_PF13354_noinserts.fasta"))
	nat = read_fasta(joinpath(folder, "PF13354_msa.fasta"))
	wnat = readdlm(joinpath(folder, "PF13354_weights.dat")) |> vec |> normalize
	R20 = read_fasta(joinpath(folder, "PSE1_rnd20_aligned_PF13354_noinserts.fasta"))
end;

# ╔═╡ f6b8e8a4-cae5-4f52-be70-ccf08c899479
function most_frequent(X)
	# most frequent element in array
	Y = countmap(X)
	return findmax(Y)
end

# ╔═╡ e9dc2b8d-b223-45af-bf36-f955f6d2dc92
begin
	tmp(T, L) = Vector{Union{Missing, T}}(undef, L)
	@kwdef mutable struct PosInfo
		pos::Int
		strats = [:wt, :cons, :iqtree, :arnet]
		stratnames = ["wild-type", "consensus", "iqtree", "autoregressive"]
		rec_int = tmp(Int, length(strats))
		rec_aa = tmp(Char, length(strats))
		fraction = tmp(Float64, length(strats))
		f_nat = tmp(Float64, length(strats))
		f_R20 = tmp(Float64, length(strats))
		f_iqtree = tmp(Float64, length(strats))
		f_arnet = tmp(Float64, length(strats))
	end
end

# ╔═╡ 348b4d28-6bd0-4191-9774-3a8660b7810d
function table_from_pi(X::PosInfo)
	L = length(X.strats)
	Y = Matrix{Any}(undef, 6, L+1)
	Y[:, 1] .= [
		"Pos. $(X.pos)", "Amino acid", "π_nat", "π_R20", "π_JTT", "π_arnet"
	]
	sigdigits=2
	for (i, s) in enumerate(X.stratnames)
		Y[:, i+1] .= [
			s, 
			string(X.rec_aa[i]), 
			"$(round(X.f_nat[i]; sigdigits))", 
			"$(round(X.f_R20[i]; sigdigits))",
			"$(round(X.f_iqtree[i]; sigdigits))", 
			"$(round(X.f_arnet[i]; sigdigits))",
		]
	end
	return Y
end

# ╔═╡ 554bd798-036c-4d09-97e2-dbffdcbd52d5
function markdown_table_from_pi(X::PosInfo)
	Y = table_from_pi(X)
	M, N = size(Y)
	S = ""
	for i in 1:M, j in 1:N
		S *= "|$(Y[i,j])"
		if j < N
			S *= ""
		else
			S *= "|\n" 
		end
		if (i == 1) && (j == N)
			for k in 1:N
				S *= "|------"
				k == N && (S*= "|\n")
			end
		end
	end
	return S
end

# ╔═╡ 6072c6f8-cd80-46ac-bf4f-ea11061be9bf
string_to_vec(s) = readdlm(IOBuffer(s[2:end-1]), ',', Int) |> vec

# ╔═╡ 613db929-1c9f-410c-8db3-5b8a72c86d46
function reconstructions_at_pos(df, strat, pos, default = wt[1][pos])
	# sub function below: 
	# if pos is in the list of positions with an error, get the corresponding rec
	# otherwise, return the default, which should be the wild-type
	function f(row)
		i = findfirst(==(pos), row.pos)
		isnothing(i) ? default : row.rec[i]
	end
	@chain df begin
		@transform @byrow :rec = string_to_vec(cols(Symbol(:rec_aa_, strat))) 
		@transform @byrow :pos = string_to_vec(cols(Symbol(:pos_, strat))) 
		map(f, eachrow(_))
	end
end

# ╔═╡ bbbabe9f-2cd0-44e3-9dd0-df6bf00387e9
function reconstructions!(df, posinfo)
	pos = posinfo.pos
	for (i, strat) in enumerate(posinfo.strats)
		if strat == :wt
			posinfo.rec_int[i] = wt[1][pos]
			posinfo.fraction[i] = 1.
		else
			frac, state = reconstructions_at_pos(df, strat, pos) |> most_frequent
			posinfo.rec_int[i] = state
			posinfo.fraction[i] = frac/100
		end
	end
	posinfo.rec_aa .= map(alphabet, posinfo.rec_int)
end

# ╔═╡ 22989241-3c55-4508-9ab0-f613c6d6687e
function get_error_positions(strat)
	X = @chain df begin
		@subset :M .== M_ref
		@select :X = cols(Symbol(:pos_, strat))
		@transform :X = map(string_to_vec, :X)
		flatten(:X)
		countmap(_.X)
		filter(x -> x[2] > 5, _)
	end
end

# ╔═╡ 21b84cb7-8e8e-4973-83e6-8e9882adfb9b
positions = let
	# finding interesting positions to look at
	err_pos = Dict(
		strat => get_error_positions(strat) for strat in strategies
	)
	interesting_positions = mapreduce(vcat, strategies) do s 
		@chain err_pos[s] keys collect
	end |> unique |> sort
end

# ╔═╡ 7d2727f0-8297-4bf4-ba25-dee4e9d6b498
for pos in positions
	@info "Position $pos - Wt $(alphabet(wt[1][pos]))"
	for strat in strategies
		recs = reconstructions_at_pos(df, strat, pos) |> countmap
		S = ""
		for (aa, n) in recs
			n > 0 && (S *= "$(alphabet(aa)): $(n/100) - ")
		end
		@info "$strat -- $S"
	end
	@info ""
end

# ╔═╡ a3bd4a20-9f10-4a39-9c08-26c3629f87c3
let
	pi = PosInfo(pos = positions[1])
	reconstructions!(df, pi)
	pi
end

# ╔═╡ 2eec84dd-803c-4957-bdeb-0ccf69a378c4
positions

# ╔═╡ a4cce3e0-bb17-4768-b10d-c052b2768a88
function parse_iqtree_model_file(model)
    lines = readlines(srcdir("iqtree_models.txt"))
    istart = findfirst(s -> occursin("model $(model)=", s), lines) + 1
    iend = findfirst(s -> occursin(";", s), lines[istart:end]) + istart - 1
    # [istart, iend-1] --> lower diag of rate matrix (not Q but Q/π)
    # [iend] --> equilibrium frequencies
    q = 20
    R = zeros(q, q)
    for (a, l) in enumerate(lines[istart:iend-1])
        R[a+1, 1:a] .= map(x -> parse(Float64, x), split(l, " "))
    end
    R .= R + R'

    p = map(x -> parse(Float64, x), split(lines[iend][1:end-1], " "))
    return R, p
end

# ╔═╡ baabeaa5-74ec-403f-a541-3230cfa4ffe5
begin
	f1_nat = site_specific_frequencies(nat, wnat; as_vec=false)
	f1_R20 = site_specific_frequencies(R20; as_vec=false)
	f1_iqtreemodel = vcat(0, parse_iqtree_model_file("JTT")[2])
end

# ╔═╡ 3fb4da85-71e1-4412-a79c-3175710c9593
position_data = map(positions) do pos
	PI = PosInfo(; pos)
	reconstructions!(df, PI)
	PI.f_nat .= f1_nat[PI.rec_int, PI.pos]
	PI.f_R20 .= f1_R20[PI.rec_int, PI.pos]
	PI.f_iqtree .= f1_iqtreemodel[PI.rec_int]
	PI.f_arnet .= ar_local_field(wt[1], PI.pos, arnet)[PI.rec_int]
	PI
end;

# ╔═╡ 59dcee32-dcbb-46ac-af89-4ec887594d06
Markdown.parse(markdown_table_from_pi(position_data[5]))

# ╔═╡ Cell order:
# ╠═d920b090-24a5-11ef-0ef8-8bd20e780e52
# ╠═4547b9d5-74b6-4d9e-8cab-6a9919ae3155
# ╠═2100eba7-5edb-41b2-be3c-9465afb3cec2
# ╠═815f6c9c-921c-40fa-bc90-4a8b26e79cdb
# ╠═fdce6423-651f-4e3b-8737-f4191439c5b1
# ╠═624397f5-9159-4482-a401-8cf71d1ddc55
# ╠═6c76669d-2c05-4c5e-bb8d-648278e91c4c
# ╠═eef3db46-fce9-45bd-a906-67fa233cdce6
# ╠═eaf5c62b-df07-45e1-8147-9581f36a66d5
# ╠═baabeaa5-74ec-403f-a541-3230cfa4ffe5
# ╟─34ee72a4-6e8f-4479-835b-c1cd6a5d6dfb
# ╠═21b84cb7-8e8e-4973-83e6-8e9882adfb9b
# ╠═7d2727f0-8297-4bf4-ba25-dee4e9d6b498
# ╠═613db929-1c9f-410c-8db3-5b8a72c86d46
# ╠═bbbabe9f-2cd0-44e3-9dd0-df6bf00387e9
# ╠═a3bd4a20-9f10-4a39-9c08-26c3629f87c3
# ╠═3fb4da85-71e1-4412-a79c-3175710c9593
# ╠═348b4d28-6bd0-4191-9774-3a8660b7810d
# ╠═554bd798-036c-4d09-97e2-dbffdcbd52d5
# ╠═2eec84dd-803c-4957-bdeb-0ccf69a378c4
# ╠═59dcee32-dcbb-46ac-af89-4ec887594d06
# ╟─fa9e3a3e-afc9-44f6-b09e-739442243d44
# ╠═9ff151f8-6807-4fef-8127-073f34a41285
# ╠═ad72228a-3381-4eef-9c15-75afa42694d0
# ╠═f6b8e8a4-cae5-4f52-be70-ccf08c899479
# ╠═e9dc2b8d-b223-45af-bf36-f955f6d2dc92
# ╠═6072c6f8-cd80-46ac-bf4f-ea11061be9bf
# ╠═22989241-3c55-4508-9ab0-f613c6d6687e
# ╠═a4cce3e0-bb17-4768-b10d-c052b2768a88
