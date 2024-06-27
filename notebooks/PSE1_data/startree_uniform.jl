### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 4f53c52c-1c0e-11ef-2eb7-77c2b7f19556
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__)

	using AutoRegressiveASR
	using BioSequenceMappings
	using Chain
	using Distributions
	using Plots
	using StatsBase
end

# ╔═╡ fd7fe7e1-25a9-4fa5-8658-895790e67791
base_folder = datadir("Stiffler/subalignments/seqstyle_uniform_treestyle_star_M1000")

# ╔═╡ b717620c-b108-4817-9f1d-11fc6453c915
rep = "id3"

# ╔═╡ c366ce00-eee4-4cac-90cb-ed97edcd1fe5
folder = joinpath(base_folder, rep)

# ╔═╡ c6767794-744f-4b34-9b8b-a89caa7507fa
readdir(folder)

# ╔═╡ a5200d7f-2787-47b5-9d89-19cd3f24c045
begin
	wt = read_fasta(joinpath(folder, "wt.fasta"))
	alphabet = Alphabet(wt)
	wt = alphabet(wt[1])
end

# ╔═╡ c3907c62-e428-4335-b3ea-d51ecb7153d9
gap_positions = findall(==('-'), wt)

# ╔═╡ a46480f1-5582-478a-b8b3-1e852b1a3fa1
begin
	rec_arnet = @chain joinpath(folder, "arnet_scale/ML/ml.fasta") begin
		read_fasta
		first
		alphabet
	end

	rec_iqtree = @chain joinpath(folder, "iqtree/IQTREE.state") begin
		ASRU.parse_iqtree_state_file
		get(_, "root", nothing)
		_.ml_seq
	end

	rec_cons = @chain joinpath(folder, "consensus/consensus.fasta") begin
		read_fasta
		first
		alphabet
	end
end

# ╔═╡ 92a96d02-7e46-4186-90f4-ad2092f2cec5
findall(i -> wt[i] != rec_arnet[i], 1:length(wt))

# ╔═╡ 45d8f7ea-3b07-43be-9d43-d0620f5fc585
findall(i -> wt[i] != rec_iqtree[i], 1:length(wt))

# ╔═╡ 6a3f2b31-49b4-4129-988f-5faa666f0174
findall(i -> wt[i] != rec_cons[i], 1:length(wt))

# ╔═╡ Cell order:
# ╠═4f53c52c-1c0e-11ef-2eb7-77c2b7f19556
# ╠═fd7fe7e1-25a9-4fa5-8658-895790e67791
# ╠═c366ce00-eee4-4cac-90cb-ed97edcd1fe5
# ╠═c6767794-744f-4b34-9b8b-a89caa7507fa
# ╠═a5200d7f-2787-47b5-9d89-19cd3f24c045
# ╠═c3907c62-e428-4335-b3ea-d51ecb7153d9
# ╠═a46480f1-5582-478a-b8b3-1e852b1a3fa1
# ╠═b717620c-b108-4817-9f1d-11fc6453c915
# ╠═92a96d02-7e46-4186-90f4-ad2092f2cec5
# ╠═45d8f7ea-3b07-43be-9d43-d0620f5fc585
# ╠═6a3f2b31-49b4-4129-988f-5faa666f0174
