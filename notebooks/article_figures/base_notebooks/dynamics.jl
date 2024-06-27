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

# ╔═╡ e53b1a24-f678-11ee-23d7-a95f62430963
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using AutoRegressiveASR
	using CSV
	using DataFrames
	using DataFramesMeta
	using PlutoUI
	using StatsBase
	using StatsPlots
end

# ╔═╡ b2b90692-ff75-41dc-aa28-6bddb70593f9
include(scriptsdir("families.jl"))

# ╔═╡ 9de44fae-5556-48a2-bc1f-fa8fa6958150
savepath = datadir("dynamics")

# ╔═╡ ca10e184-17ac-4f41-9abc-3ac104325c19
MOD = pluto_ingredients(scriptsdir("misc/hamming_v_time_potts.jl"))

# ╔═╡ 6778b945-efe1-47fc-8524-e49447a1aed1
fam_picker = @bind family Select(
	[v => k for (k, v) in families]
)

# ╔═╡ 7d038ac0-4c3e-4d4c-a78f-33ad8bdeda14
data_potts, data_ardca = let
	dat, _ = produce_or_load(
		family, savepath;
		filename = X -> X["prefix"], 
		prefix = "hamming_vs_time_potts_arnet", 
		verbose=true,
	) do fam
        data_potts = MOD.hamming_v_time_potts(fam; Nsamples = 500)
        data_ardca = MOD.hamming_v_time_ardca(fam; Nsamples = 500)
        Dict(
            "family" => fam,
            "data_potts" => data_potts,
            "data_ardca" => data_ardca,
            "timestamp" => now()
        )
	end
	dat["data_potts"], dat["data_ardca"]
end;

# ╔═╡ 6a4f00b5-d5ac-4634-b088-681b1f9cdab0
md"# Utils"

# ╔═╡ 92e6a51a-c3e5-4a59-a1bc-2ecba3a17d4f
function geometric_middle(x1, x2)
	if x1 > x2
		return geometric_middle(x2, x1)
	elseif x1 == 0 
		return (x2-x1)/2
	elseif x1 == x2
		return x1
	end

	α = x2/x1
	return x1 * sqrt(α)
end

# ╔═╡ 376e70f4-764a-4f1c-8231-69db94bf26b7
begin
	function interp_hamming(htarget, hvals, tvals)
		ilow = findlast(<=(htarget), hvals)
		ihigh = findfirst(>(htarget), hvals)
		isnothing(ihigh) && return tvals[ilow]
		
		α = (htarget - hvals[ilow]) / (hvals[ihigh] - hvals[ilow])
		return (1-α)*tvals[ilow] + α*tvals[ihigh]
	end
	function interp_hamming(htarget, data::DataFrame) 
		return @df data interp_hamming(htarget, :hamming_to_init_av, :t)
	end
end

# ╔═╡ 372373fd-6800-4781-be8b-2208baddcb42
hvals, time_potts, time_ardca = let
	hmax = 0.95 * min(
		maximum(data_potts.hamming_to_init_av), 
		maximum(data_ardca.hamming_to_init_av)
	)
	hvals = collect(range(0, hmax, 100))

	(
		hvals, 
		map(h -> interp_hamming(h, data_potts), hvals), 
		map(h -> interp_hamming(h, data_ardca), hvals),
	)
end

# ╔═╡ 00e601d7-357f-44e1-87a7-53d0dd9cf162
let p = plot()
	@df data_potts plot!(
		:t, :hamming_to_init_av; 
		ribbon = :hamming_to_init_conf95, label = "",
	)
	plot!(
		xlim = extrema(time_potts[2:end]),
		xscale = :log10,
	)
end

# ╔═╡ 18294804-8bb0-4b2a-ac88-b992cb7775fb
let p = plot()
	@df data_ardca plot!(
		:t, :hamming_to_init_av; 
		ribbon = :hamming_to_init_conf95, label = "",
	)
	plot!(
		xlim = extrema(time_ardca[2:end]),
		xscale = :log10,
	)
end

# ╔═╡ 15189e13-4dcd-4162-8416-00e547bef791
hvals

# ╔═╡ d84760d7-c5e1-4d42-8df5-eb7b342890d7
let p = plot()
	plot(time_ardca, time_potts, label="")

	T = 1.25
	i = findmin(t -> abs(t-T), time_ardca)[2]
	T_potts = round(time_potts[i]; sigdigits=3)
	
	plot!([0, time_ardca[i]], [T_potts, T_potts], line=(:black, :dash), label="")
	plot!([T, T], [0, time_potts[i]], line=(:black, :dash), label="")
	annotate!(0, T_potts * .9, text("$(T_potts) sweeps", :left))

	plot!(
		title = family["prefix"],
		xlabel = "time - ArDCA",
		ylabel = "sweeps - Potts",
	)
end

# ╔═╡ 0ddb2e7b-1c3a-4577-a013-68308b1d10fc
let p1 = plot(), p2 = plot()
	plot!(p1, time_ardca, hvals)
	@df data_ardca scatter!(p1, :t, :hamming_to_init_av)
	
	plot!(p2, time_potts, hvals)
	@df data_potts scatter!(p2, :t, :hamming_to_init_av)

	plot(p1, p2, layout=grid(2,1), size = (400, 600))
	plot!(title = "Should fit - Check that interp is working well")
end

# ╔═╡ b1b588ed-6e90-4b97-9708-dba94433007f
interp_hamming(.5, data_ardca)

# ╔═╡ f476a055-3f0d-4230-9b47-31ef17a81f4a
let
	@df data_ardca plot(:t, :hamming_to_init_av, marker=:o)
	@df data_potts plot!(:t / 200, :hamming_to_init_av, marker=:o)
	h = .42
	hline!([h])
	t = interp_hamming(h, data_ardca)
	vline!([t])
	
	plot!(
		xlim = (0, 3),
		# xscale=:log10,
	)
end

# ╔═╡ Cell order:
# ╠═e53b1a24-f678-11ee-23d7-a95f62430963
# ╠═b2b90692-ff75-41dc-aa28-6bddb70593f9
# ╠═9de44fae-5556-48a2-bc1f-fa8fa6958150
# ╠═ca10e184-17ac-4f41-9abc-3ac104325c19
# ╠═7d038ac0-4c3e-4d4c-a78f-33ad8bdeda14
# ╠═372373fd-6800-4781-be8b-2208baddcb42
# ╠═00e601d7-357f-44e1-87a7-53d0dd9cf162
# ╠═18294804-8bb0-4b2a-ac88-b992cb7775fb
# ╠═6778b945-efe1-47fc-8524-e49447a1aed1
# ╠═15189e13-4dcd-4162-8416-00e547bef791
# ╠═b1b588ed-6e90-4b97-9708-dba94433007f
# ╠═d84760d7-c5e1-4d42-8df5-eb7b342890d7
# ╟─0ddb2e7b-1c3a-4577-a013-68308b1d10fc
# ╟─6a4f00b5-d5ac-4634-b088-681b1f9cdab0
# ╠═92e6a51a-c3e5-4a59-a1bc-2ecba3a17d4f
# ╠═376e70f4-764a-4f1c-8231-69db94bf26b7
# ╠═f476a055-3f0d-4230-9b47-31ef17a81f4a
