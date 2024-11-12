### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 38c89bb4-986b-11ef-36c3-35a8cf357289
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__)
end

# ╔═╡ e6b895da-cb04-4242-865b-fd7ac3e09c4a
begin
	using AutoRegressiveASR
	using AncestralSequenceReconstruction
	using Chain
	using FASTX
	using JLD2
	using Plots
	using StatsBase
	using Statistics
	using TreeTools
end

# ╔═╡ 9c4c2b3e-4843-4582-ae2a-e9fa754d3a0a
include(scriptsdir("families.jl"))

# ╔═╡ 14143d83-9270-4588-8c6e-fd76df94e5ba
md"""
# Reading data
"""

# ╔═╡ a17cb240-d40e-4914-85a6-58929edc473d
maindir = datadir(
	"simulated/_tests",
	"PF00072_nleaves=100_ntrees=1_opt_bl=fromreal_treeheight=2.0",
	"data/1",
)

# ╔═╡ a6955baa-790e-4c7e-aa26-110fa88b888b
begin
	newick_real = joinpath(maindir, "tree.nwk")
	fasta_simulated = joinpath(maindir, "alignment_leaves.fasta")
	fasta_nat = families["PF00072"]["aln_nat_small"]
end

# ╔═╡ 84afa2e8-a800-4020-bd73-90566abb4f00
tree_real = let # Construct real tree with sequences mapped to leaves
    seqmap = FASTAReader(open(fasta_simulated, "r")) do reader
        map(rec -> identifier(rec) => sequence(rec),reader)
    end
	L = length(first(seqmap)[2])
	q = 21
	T() = ASR.AState{q}(;L)
	tree = read_tree(newick_real; node_data_type = T)
	ASR.sequences_to_tree!(tree, seqmap; alphabet=ASR.Alphabet(:aa))
	tree
end

# ╔═╡ 9bbdd0ab-4612-4726-b1e8-32ad8398edc9
begin
	@load families["PF00072"]["arnet"] arnet
	global_profile = ASR.ProfileModel(fasta_nat; pc=1e-2)
	local_profile = ASR.ProfileModel(
		fasta_simulated; 
		pc=0.1, reweighting=true, θ=0.2,
	)
end

# ╔═╡ 298e758d-40cc-4e1f-9f0b-9443312d8128
md"""
# Optim
"""

# ╔═╡ 161fe486-d3e7-4f13-9ac9-d6f2927fed7e
rconv = 1e-2

# ╔═╡ 3f17bf98-c1a1-42e1-a469-866dc9256183
opt_strat = ASRMethod(; 
	joint=false, optimize_branch_length_cycles=2, verbosity=2,
);

# ╔═╡ 6d47c41c-66e4-4cdf-8c2d-b5b1e4b1d984
tree_globalprofile, lk_globalprofile = ASR.optimize_branch_length(
	tree_real, global_profile, opt_strat; rconv
)

# ╔═╡ d4493f7d-8fb9-4484-9d00-4b075b54ce80
tree_localprofile, lk_localprofile = ASR.optimize_branch_length(
	tree_real, local_profile, opt_strat; rconv
)

# ╔═╡ ffde2629-8367-41d2-ad12-787edd5c6c94
let p = plot()
	dlocal = map(branch_length, tree_localprofile) |> skipmissing
	dglobal = map(branch_length, tree_globalprofile) |> skipmissing
	scatter!(dglobal, dlocal)
	plot!([0, maximum(dglobal)], [0, maximum(dglobal)])
	# maximum(dglobal)
end

# ╔═╡ 6d1261c1-aea1-4c98-a222-56dc51e01ff8


# ╔═╡ a9e298ff-b20f-409d-820b-041aa8a14f6a
map(label, children(tree_localprofile.root))

# ╔═╡ 39144c92-0a0e-4163-9c22-ab6a5201b478
@chain tree_localprofile begin
	nodes
	collect
	nds = filter(x -> !isroot(x), _)
	findall(x -> branch_length(x) > 2, _)
	nds[_]
	map(label, _)
end

# ╔═╡ 5414e3a2-35cc-422c-8be1-d4807c0259a6
clade_seqs, nonclade_seqs = let
	node = tree_real["internal_95"]
	clade = collect(POTleaves(node))
	nonclade = filter(!in(clade), collect(leaves(tree_real)))
	
	(
		mapreduce(x -> data(x).sequence, hcat, clade),
		mapreduce(x -> data(x).sequence, hcat, nonclade),
	)
end;

# ╔═╡ 35127be5-1590-4cc8-abfd-f15b15e5552e
i=28

# ╔═╡ e4ef9069-2432-483b-8cc2-d3ff055054f4
countmap(clade_seqs[i, :])

# ╔═╡ 13c05ad7-c049-4b8e-9483-a4005bc3fa8b
countmap(nonclade_seqs[i, :])

# ╔═╡ 4d3a4e34-923b-48a3-9383-b6e0a92b02e1
local_profile.P[1]

# ╔═╡ bd93638d-5859-455d-bccd-97f63459ada9
md"""
# Utils
"""

# ╔═╡ 341e031f-c463-44a9-a557-52e92b4d0ace
function plt_distance_matrix(tree_real, trees...)
	p = plot()
	Dreal = TreeTools.distance_matrix(tree_real) |> vec
	for (name, tree) in trees
	   D = TreeTools.distance_matrix(tree) |> vec
	   scatter!(Dreal, D, label=name, marker=(5, stroke(0)))
	end
	plot!(
	   [0,maximum(Dreal)], [0, maximum(Dreal)]; 
	   line=(:black, :dash), label="",
	)
	plot!(
		xlabel="Real", 
		ylabel="Inferred",
	)
end

# ╔═╡ 90a2c345-0e29-42ee-95c4-a51204c85d4c
plt_distance_matrix(
	tree_real, 
	("global profile", tree_globalprofile),
	("local profile", tree_localprofile),
)

# ╔═╡ dcd0d656-ae13-4042-89a0-3daeac60af82
function plt_branch_length(tree_real, trees...)
	p = plot()
	nodes = filter(!=("root"), map(label, POT(tree_real)))
	
	Dreal = map(x -> branch_length(tree_real[x]), nodes)
	for (name, tree) in trees
	   D = map(x -> branch_length(tree[x]), nodes)
	   scatter!(Dreal, D, label=name, marker=(5, stroke(0)))
	end
	plot!(
	   [0,maximum(Dreal)], [0, maximum(Dreal)]; 
	   line=(:black, :dash), label="",
	)
	plot!(
		xlabel="Real", 
		ylabel="Inferred",
	)
end

# ╔═╡ da60bc02-ba03-468c-aa1b-6bd2ee7ef533
plt_branch_length(
	tree_real, 
	("global profile", tree_globalprofile),
	("local profile", tree_localprofile),
)

# ╔═╡ Cell order:
# ╠═38c89bb4-986b-11ef-36c3-35a8cf357289
# ╠═e6b895da-cb04-4242-865b-fd7ac3e09c4a
# ╟─14143d83-9270-4588-8c6e-fd76df94e5ba
# ╠═9c4c2b3e-4843-4582-ae2a-e9fa754d3a0a
# ╠═a17cb240-d40e-4914-85a6-58929edc473d
# ╠═a6955baa-790e-4c7e-aa26-110fa88b888b
# ╠═84afa2e8-a800-4020-bd73-90566abb4f00
# ╠═9bbdd0ab-4612-4726-b1e8-32ad8398edc9
# ╟─298e758d-40cc-4e1f-9f0b-9443312d8128
# ╠═161fe486-d3e7-4f13-9ac9-d6f2927fed7e
# ╠═3f17bf98-c1a1-42e1-a469-866dc9256183
# ╠═6d47c41c-66e4-4cdf-8c2d-b5b1e4b1d984
# ╠═d4493f7d-8fb9-4484-9d00-4b075b54ce80
# ╠═90a2c345-0e29-42ee-95c4-a51204c85d4c
# ╠═da60bc02-ba03-468c-aa1b-6bd2ee7ef533
# ╠═ffde2629-8367-41d2-ad12-787edd5c6c94
# ╠═6d1261c1-aea1-4c98-a222-56dc51e01ff8
# ╠═a9e298ff-b20f-409d-820b-041aa8a14f6a
# ╠═39144c92-0a0e-4163-9c22-ab6a5201b478
# ╠═5414e3a2-35cc-422c-8be1-d4807c0259a6
# ╠═35127be5-1590-4cc8-abfd-f15b15e5552e
# ╠═e4ef9069-2432-483b-8cc2-d3ff055054f4
# ╠═13c05ad7-c049-4b8e-9483-a4005bc3fa8b
# ╠═4d3a4e34-923b-48a3-9383-b6e0a92b02e1
# ╟─bd93638d-5859-455d-bccd-97f63459ada9
# ╠═341e031f-c463-44a9-a557-52e92b4d0ace
# ╠═dcd0d656-ae13-4042-89a0-3daeac60af82
