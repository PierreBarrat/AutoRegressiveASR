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

# ╔═╡ 157c8ab0-d63d-11ee-3794-27a66b691088
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using AutoRegressiveASR
	using CSV
	using DataFrames
	using DataFramesMeta
	using DCATools
	using JSON3
	using PlutoUI
	using StatsBase
	using StatsPlots
end

# ╔═╡ 7ebf97e2-fcca-4420-a01b-835de467ef89
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ c6d4a813-23f9-48a3-8418-398d94886da4
begin
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 0ab653c4-95e2-4676-bfda-d5b15f3d500f
_fs = @bind folder_full Select(readdir(datadir("simulated/potts_yule"); join=true))

# ╔═╡ f9b87a61-f1cd-4e57-a6c3-d2e6733a12ff
folder = basename(folder_full)

# ╔═╡ 6ebfeb60-837f-45e7-a918-bcbc91add5a2
title_id = split(folder, "_")[1]

# ╔═╡ 692ed38d-95b0-4003-b337-27c62efc9973
md"## Loading data"

# ╔═╡ 577be925-516c-4efb-bbe0-3a07d6dfd2e0
simulation_parameters = JSON3.read(
	joinpath(folder_full, "simulation_parameters.json"), Dict
)

# ╔═╡ 39a98d0e-d550-4329-8048-b0462c8bd2e2
generative_model = simulation_parameters["potts_file"] |> DCAGraph

# ╔═╡ f4c24956-9f41-4b45-9617-ca965bcc7aab
model_consensus = let
	sample_file = projectdir(simulation_parameters["sample_potts_file"])
	Seq = if isnothing(sample_file) || !isfile(sample_file)
		T = 100
		M = 1000
		@info "sampling generative model for consensus. Check eq. time (default $T)"
		S = DCATools.sample(generative_model, M; Twait = T)
	else
		read_msa(sample_file)
	end
	cons = DCATools.consensus(Seq) # a DCASample object
	DCATools.num_to_aa(cons[1], cons.mapping) # a string sequence
end

# ╔═╡ e0922d49-043e-494c-947e-e58e3a808224
md"# Figures"

# ╔═╡ 63b87289-1634-44ac-8e2c-e51ac3a7c271
md"## Hamming distance with gaps"

# ╔═╡ 471f444a-7de5-4399-baba-8d189ed06670
begin
	# smoothing width
	w = 15
	smoothing_alg = :hist
	pal = palette(:default)
	strat_clr = Dict("iqtree" => pal[1], "autoregressive" => pal[2], "real" => pal[3])
end

# ╔═╡ 04ab3d51-0c30-4bd9-b08a-0c7e0540f3b1
md"## Hamming distance without gaps"

# ╔═╡ b96f65d5-cb63-4396-967a-f190e1979dcf
md"## Likelihood of reconstructed sequences"

# ╔═╡ 70b1c534-5f41-47f0-8fec-27201b0325d9
md"## Hamming distance to consensus"

# ╔═╡ 26f9e913-0b55-4d2f-9109-492787d386c6
md"## Functions"

# ╔═╡ 3fa3dd10-edbb-4ffe-94ec-55f665ba62b4
begin
	bayesian(strategies) = filter(x -> length(x)>1 && x[2]=="Bayes", strategies)
	ml(strategies) = filter(strategies) do x 
		length(x) < 2 && return false
		x[2] == "ML" || x[2] == "ml"
	end
	real(strategies) = filter(==(["real"]), strategies)
	reconstruction(strategies) = filter(!=(["real"]), strategies)
	strat_label(strat) = joinpath(strat...)

	iqtree(strategies) = filter(x -> x[1]=="iqtree", strategies)
	ar(strategies) = filter(x -> x[1]=="autoregressive", strategies)

	function label2(strat)
		length(strat) == 1 && return strat[1]
		strat[2] == "Bayes" ? "" : strat[1]
	end
	function linestyle(strat)
		lw = 4
		return if length(strat) > 1 && strat[2] == "Bayes"
			(lw, :dash, strat_clr[strat[1]])
		else
			(lw, strat_clr[strat[1]])
		end
	end
end

# ╔═╡ 9922c4b4-04ff-4bf8-8e56-8959a49e1c91
!=(2)

# ╔═╡ d2100e2a-8a86-4aad-bff0-c1420aabdac7
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ ecdbbfcf-0d17-49de-9f4f-18ae0d2e85eb
begin
	LM = ingredients(scriptsdir("figures_and_results/analyze_results_and_write_df.jl"))
	md"""
	This puts in scope the function `analyze_results_and_write(folder; strategies)`. It searches `folder` for `strat[1]/strat[2]` where `strat` is by default in 
	- `["iqtree", "ML"]
	- ...
	"""
end

# ╔═╡ a6cce528-8e5e-45b1-8019-8ac1f4f6f633
analyze_log = if isfile(joinpath(folder_full, "analyze_results_log.json"))
	JSON3.read(joinpath(folder_full, "analyze_results_log.json"), Dict)
else
	LM.analyze_results_and_write(folder_full)
	JSON3.read(joinpath(folder_full, "analyze_results_log.json"), Dict)
end;

# ╔═╡ 8cbaee45-3842-44e6-a54d-9977f734462d
data = let
	itr = zip(analyze_log["strategies"], analyze_log["filenames"])
	dat = Dict()
	for (strat, df) in itr
		file = joinpath(folder_full, df)
		dat[strat] = @subset DataFrame(CSV.File(file)) :label .!= "root"
	end
	dat
end;

# ╔═╡ 39c9733b-6bbe-4d5d-bd70-3c90d5844b13
strategies = collect(keys(data))

# ╔═╡ b93a2a13-ace0-4d1c-9134-57d9e9ff0207
map(label2, strategies)

# ╔═╡ 594bf078-4cf3-404a-bc55-d208fc4d8394
strat_label(strategies[1])

# ╔═╡ 6d13312e-f1cf-4d25-9eaa-6784ee6386f6
let p = plot()
	for (i, strat) in enumerate(bayesian(strategies))
		# @df data[strat] scatter!(
		# 	:node_depth, :hamming_to_real;
		# 	label=strat, color = i, marker=(3, .25, stroke(0))
		# )
		x, y = ASRU.easy_smooth(
            data[strat], :node_depth, :hamming_to_real; w, alg=smoothing_alg,
        )
		plot!(x, y, label=strat_label(strat), color=i, line=(3))
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real - Bayes",
		frame = :box,
		legend = :bottomright,
	)
end

# ╔═╡ 71613c45-8a7d-42ea-88bc-b4d5132f6a65
let p = plot()
	for (i, strat) in enumerate(ml(strategies))
		x, y = ASRU.easy_smooth(
            data[strat], :node_depth, :hamming_to_real; w, alg=smoothing_alg,
        )
		plot!(x, y, label=strat_label(strat), color=i, line=(3))
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real - ML",
		frame = :box,
		legend = :bottomright,
	)
end

# ╔═╡ 103932ee-7c76-44f0-a5cb-78a141687996
plt_hamming_to_real_gaps = let p = plot()
	for strat in reconstruction(strategies)
		# @df data[strat] scatter!(
		# 	:node_depth, :hamming_to_real;
		# 	label="", color = strat_clr[strat[1]], marker=(3, .25, stroke(0))
		# )
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_real; w, alg=smoothing_alg,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real",
		frame = :box,
		legend = :bottomright,
	)
end

# ╔═╡ 8f6a4029-04e8-40f8-8dcd-144d5ede6c77
plt_hamming_to_real_nogaps = let p = plot()
for strat in reconstruction(strategies)
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_real_nogap; w, alg=smoothing_alg,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Hamming distance to real",
		title = "$title_id - Hamming to real without gaps",
		frame = :box,
		legend = :bottomright,
	)
end

# ╔═╡ d13b857a-4071-4f96-816e-ea7949b3581c
plt_likelihood = let p = plot()
	for strat in strategies
		# @df data[strat] scatter!(
		# 	:node_depth, :loglikelihood;
		# 	label="", color = strat_clr[strat[1]], marker=(3, .25, stroke(0))
		# )
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :loglikelihood; w=25, alg=:hist,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "Likelihood (teacher model)",
		title = "$title_id - Likelihood",
		frame = :box,
		legend = false#:topleft,
	)
	p
end

# ╔═╡ 41dbf0f5-8471-4d7d-901d-dd63b63783d2
plt_hamming_to_aln_consensus = let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_aln_consensus; w, alg=smoothing_alg,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end

	# 
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "",
		title = "Hamming to alignment consensus",
		frame = :box,
		legend = :bottomleft,
	)
end

# ╔═╡ f0a26907-0304-4841-9401-9628a41a7c72
plt_hamming_to_model_consensus = let p = plot()
	for strat in strategies
		x, y = ASRU.easy_smooth(
			data[strat], :node_depth, :hamming_to_model_consensus; w, alg=smoothing_alg,
		)
		plot!(x, y, label=label2(strat), line=linestyle(strat))
	end

	# 
	
	plot!(
		xlabel = "Node depth (mcmc steps)",
		ylabel = "",
		title = "Hamming to generative model consensus",
		frame = :box,
		legend = :bottomleft,
	)
end

# ╔═╡ Cell order:
# ╠═157c8ab0-d63d-11ee-3794-27a66b691088
# ╠═7ebf97e2-fcca-4420-a01b-835de467ef89
# ╠═c6d4a813-23f9-48a3-8418-398d94886da4
# ╠═0ab653c4-95e2-4676-bfda-d5b15f3d500f
# ╠═f9b87a61-f1cd-4e57-a6c3-d2e6733a12ff
# ╠═6ebfeb60-837f-45e7-a918-bcbc91add5a2
# ╟─692ed38d-95b0-4003-b337-27c62efc9973
# ╠═577be925-516c-4efb-bbe0-3a07d6dfd2e0
# ╠═39a98d0e-d550-4329-8048-b0462c8bd2e2
# ╠═f4c24956-9f41-4b45-9617-ca965bcc7aab
# ╠═a6cce528-8e5e-45b1-8019-8ac1f4f6f633
# ╠═8cbaee45-3842-44e6-a54d-9977f734462d
# ╠═39c9733b-6bbe-4d5d-bd70-3c90d5844b13
# ╟─e0922d49-043e-494c-947e-e58e3a808224
# ╟─63b87289-1634-44ac-8e2c-e51ac3a7c271
# ╟─6d13312e-f1cf-4d25-9eaa-6784ee6386f6
# ╟─71613c45-8a7d-42ea-88bc-b4d5132f6a65
# ╠═471f444a-7de5-4399-baba-8d189ed06670
# ╠═103932ee-7c76-44f0-a5cb-78a141687996
# ╟─04ab3d51-0c30-4bd9-b08a-0c7e0540f3b1
# ╟─8f6a4029-04e8-40f8-8dcd-144d5ede6c77
# ╟─b96f65d5-cb63-4396-967a-f190e1979dcf
# ╟─d13b857a-4071-4f96-816e-ea7949b3581c
# ╟─70b1c534-5f41-47f0-8fec-27201b0325d9
# ╟─41dbf0f5-8471-4d7d-901d-dd63b63783d2
# ╟─f0a26907-0304-4841-9401-9628a41a7c72
# ╟─26f9e913-0b55-4d2f-9109-492787d386c6
# ╟─ecdbbfcf-0d17-49de-9f4f-18ae0d2e85eb
# ╠═3fa3dd10-edbb-4ffe-94ec-55f665ba62b4
# ╠═b93a2a13-ace0-4d1c-9134-57d9e9ff0207
# ╠═594bf078-4cf3-404a-bc55-d208fc4d8394
# ╠═9922c4b4-04ff-4bf8-8e56-8959a49e1c91
# ╠═d2100e2a-8a86-4aad-bff0-c1420aabdac7
