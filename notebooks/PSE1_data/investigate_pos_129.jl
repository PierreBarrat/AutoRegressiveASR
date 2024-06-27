### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 910fdb5a-24dc-11ef-3c86-af9f3dd2e047
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__)

	using ArDCA
	using AncestralSequenceReconstruction
	using BioSequenceMappings
	using Chain
	using CSV
	using DataFrames
	using DataFramesMeta
	using DelimitedFiles
	using IterTools
	using JLD2
	using StatsPlots
	using StatsBase
	using TreeTools
end

# ╔═╡ 5e7ac782-9306-4597-9cf1-b74a75084f2b
folder = datadir(
	"Stiffler/subalignments/seqstyle_uniform_treestyle_star_nobinarize/M160/id1"
)

# ╔═╡ 86e2faf5-bf64-4092-95ab-dc3f33a4cc16
pos = 129

# ╔═╡ bd5d3c9f-7803-405c-a129-56224ca13ef9
alphabet = BioSequenceMappings.Alphabet(:aa)

# ╔═╡ 9385ca4e-a142-4856-8b2f-2eeb5b7da8ba
begin
	arnet = JLD2.load(
		datadir("Stiffler/arnet/arnet_PF13354_lJ0.01_lH0.001.jld2")
	)["arnet"]
	evo_arnet = ASR.AutoRegressiveModel(arnet)
end

# ╔═╡ d18289a4-8a58-4ac1-b061-82303fc0c187
strategy = ASRMethod(; 
	joint=false, ML=true, optimize_branch_length=false, verbosity=2
)

# ╔═╡ 6a9e8267-5838-4746-9835-875ccb0dfb6e
begin
	tree = joinpath(folder, "arnet_scale/tree_opt.nwk")
	aln_file = joinpath(folder, "R20.fasta")
end

# ╔═╡ 553422af-5b5c-48ef-a11b-60da08c23fec
f_aln = let
	aln = read_fasta(aln_file)
	f = site_specific_frequencies(aln; as_vec=false)
	f[:, pos]
end

# ╔═╡ e2a9005d-0558-47d1-947b-7bf687b23341
let
	X = ArDCA.sample(arnet, 1000)
	site_specific_frequencies(Alignment(X, alphabet); as_vec=false)[:, pos]
end

# ╔═╡ 0730becd-63dc-4959-a8a7-1ce03ea0c3e7
arnet.idxperm

# ╔═╡ d244c32b-01d5-426f-a827-94105e8a6338
out = infer_ancestral(
	tree, aln_file, evo_arnet, strategy;
)

# ╔═╡ 875212fd-7b75-4c8b-bba3-a1bd869c47f4
f_aln

# ╔═╡ 31402122-b62e-491f-8fa5-7230e36cfb08
let
	n = first(leaves(out[1]))
	n.data.pstates[pos].weights
end

# ╔═╡ b638e4f2-7137-4915-a67e-51879397fc34
map(leaves(out[1])) do n 
	n.data.pstates[pos].weights.π[11]
end |> sort

# ╔═╡ 6d5f98a5-f301-479c-a2f5-c19df0814730


# ╔═╡ 987e7d20-ac27-4a99-9806-7534ed63d680
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

# ╔═╡ 22a862ef-890d-44fd-a970-33162ff2e98f
local_field_aln = let
	aln = read_fasta(aln_file)
	mapreduce(hcat, eachsequence(aln)) do ref
		ar_local_field(ref, pos, arnet)
	end
end

# ╔═╡ 8cb4bb3b-47a7-4404-830d-4259de19a0d1
let
	local_field_aln_P = [X[alphabet('P')] for X in eachcol(local_field_aln)]
	local_field_aln_L = [X[alphabet('L')] for X in eachcol(local_field_aln)]
	@chain local_field_aln_L sort# log.(_)
end

# ╔═╡ Cell order:
# ╠═910fdb5a-24dc-11ef-3c86-af9f3dd2e047
# ╠═5e7ac782-9306-4597-9cf1-b74a75084f2b
# ╠═86e2faf5-bf64-4092-95ab-dc3f33a4cc16
# ╠═bd5d3c9f-7803-405c-a129-56224ca13ef9
# ╠═9385ca4e-a142-4856-8b2f-2eeb5b7da8ba
# ╠═d18289a4-8a58-4ac1-b061-82303fc0c187
# ╠═6a9e8267-5838-4746-9835-875ccb0dfb6e
# ╠═553422af-5b5c-48ef-a11b-60da08c23fec
# ╠═e2a9005d-0558-47d1-947b-7bf687b23341
# ╠═0730becd-63dc-4959-a8a7-1ce03ea0c3e7
# ╠═d244c32b-01d5-426f-a827-94105e8a6338
# ╠═875212fd-7b75-4c8b-bba3-a1bd869c47f4
# ╠═31402122-b62e-491f-8fa5-7230e36cfb08
# ╠═b638e4f2-7137-4915-a67e-51879397fc34
# ╠═6d5f98a5-f301-479c-a2f5-c19df0814730
# ╠═22a862ef-890d-44fd-a970-33162ff2e98f
# ╠═8cb4bb3b-47a7-4404-830d-4259de19a0d1
# ╠═987e7d20-ac27-4a99-9806-7534ed63d680
