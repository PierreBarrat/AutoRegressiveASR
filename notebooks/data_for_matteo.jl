### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 4bdf4d6c-d254-4073-95c8-1053deace0b4
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__) 
	
	using ArDCA
	using AutoRegressiveASR
	using BioSequenceMappings
	using CSV
	using DataFrames 
	using DataFramesMeta 
	using JLD2
end

# ╔═╡ 5ba5ccf5-fa6c-480f-a9b8-97c9942ac4e1


# ╔═╡ 818fec0a-2e2e-11ef-2cd2-1391b2db5606
data_path = "/home/pierrebc/Documents/BaleLabo/Data/AutoRegressiveASR/data/Stiffler/subalignments/Results"

# ╔═╡ cc41b415-969b-4c27-9473-63e189c2fc8e
df = JLD2.load(joinpath(data_path, "data.jld2"))["data_wcode"];

# ╔═╡ 8c3b0a9e-9ab4-4765-857f-68665c5d3406
md"## iqtree model"

# ╔═╡ 2bea4eea-9642-467a-90f9-117da9912568
alphabet_permutation = let
	map(enumerate(Alphabet(:aa).string)) do (i,a)
		if in(a, AutoRegressiveASR.alphabet_iqtree)
			AutoRegressiveASR.alphabet_iqtree(a)
		else
			missing
		end
	end |> skipmissing |> collect
end

# ╔═╡ 60b7766a-16a5-494e-9301-886dff3a1b5d
D_iqtree, P_iqtree = let
	D, P = AutoRegressiveASR.parse_iqtree_model_file("JTT")
	P = vcat(0, P[alphabet_permutation])
	Base.permutecols!!(D, copy(alphabet_permutation))
	Base.permutecols!!(D', copy(alphabet_permutation))'
	D, P
end

# ╔═╡ 1c5017fb-9721-4529-90d2-c8c7d95522cd
md"# Profiles in wt context"

# ╔═╡ f2498ced-c82a-4567-9bcb-1d96ca05e00c
arnet = JLD2.load(datadir("Stiffler/arnet/arnet_PF13354_lJ0.01_lH0.001.jld2"))["arnet"];

# ╔═╡ 878a10a5-5c44-4b27-ba67-60dbcaba64b4
wt = read_fasta(
	datadir("Stiffler/aligned_data_ref/PSE1_aligned_PF13354_noinserts.fasta")
)[1] |> collect

# ╔═╡ 9454c379-b61f-42be-8eb5-cd91049a4919


# ╔═╡ fb570faf-8516-4c3e-92dd-ad8d64738f5e
positions = [64, 68, 129, 147, 152, 212]

# ╔═╡ fa3e9998-b7f5-4015-a55d-e9829c14495b
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

# ╔═╡ 2a30426d-455b-4197-baa7-7b1f87907aa9
arnet_wt_profiles = Dict(i => ar_local_field(wt, i, arnet) for i in positions)

# ╔═╡ f15e851e-5283-455e-b9af-7ea1c2659d6f
md"## Save"

# ╔═╡ 147feed0-0a8b-45b5-ba35-5f506101daed
@save(
	datadir("Stiffler/subalignments/Results/iqtree_arnet_profiles.jld2"),
	alphabet = Alphabet(:aa).string, P_iqtree, D_iqtree, arnet_wt_profiles,
)

# ╔═╡ 14c7b1d6-cc3d-428a-bdfb-815f55e1ef80


# ╔═╡ Cell order:
# ╠═5ba5ccf5-fa6c-480f-a9b8-97c9942ac4e1
# ╠═4bdf4d6c-d254-4073-95c8-1053deace0b4
# ╠═818fec0a-2e2e-11ef-2cd2-1391b2db5606
# ╠═cc41b415-969b-4c27-9473-63e189c2fc8e
# ╟─8c3b0a9e-9ab4-4765-857f-68665c5d3406
# ╠═2bea4eea-9642-467a-90f9-117da9912568
# ╠═60b7766a-16a5-494e-9301-886dff3a1b5d
# ╟─1c5017fb-9721-4529-90d2-c8c7d95522cd
# ╠═f2498ced-c82a-4567-9bcb-1d96ca05e00c
# ╠═878a10a5-5c44-4b27-ba67-60dbcaba64b4
# ╠═9454c379-b61f-42be-8eb5-cd91049a4919
# ╠═fb570faf-8516-4c3e-92dd-ad8d64738f5e
# ╠═fa3e9998-b7f5-4015-a55d-e9829c14495b
# ╠═2a30426d-455b-4197-baa7-7b1f87907aa9
# ╠═f15e851e-5283-455e-b9af-7ea1c2659d6f
# ╠═147feed0-0a8b-45b5-ba35-5f506101daed
# ╠═14c7b1d6-cc3d-428a-bdfb-815f55e1ef80
