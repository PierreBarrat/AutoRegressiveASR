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

# ╔═╡ ec42999a-edaf-11ee-21a4-8558cad7bc28
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using ArDCA
	using AutoRegressiveASR
	using Chain
	using DCATools
	using EasyFit
	using JLD2
	using Measures
	using Plots
	using PlutoUI
	using StatsBase
end

# ╔═╡ 5dcb047d-46b6-45d7-9780-a98c0055eb7a
families = let
	FAM = pluto_ingredients(scriptsdir("families.jl"))
	FAM.families
end

# ╔═╡ 7df58d8b-95a4-4a2c-ad5f-a1747797598c
fam_picker = @bind family Select(map(reverse, collect(pairs(families))))

# ╔═╡ dbdf3b60-8579-472d-8236-7005788d3ceb
family["arnet"]

# ╔═╡ 559526a0-8933-4508-81e4-6a61d6e1d229
arnet = JLD2.load(family["arnet"])["arnet"]

# ╔═╡ eadcdce5-dd16-409e-8e7b-d70edc2eb380
L = length(arnet.J) + 1

# ╔═╡ 57d52ec7-6085-4c44-aeef-0bd48baa156a
q = length(arnet.p0)

# ╔═╡ 36e7cb69-22ac-4c09-bdf3-37c217b03d29
potts = DCAGraph(family["potts"])

# ╔═╡ 92b6f8a2-b01a-4a9e-9167-4a4a80b14966
aln = DCASample(ArDCA.sample(arnet, 10_000)');

# ╔═╡ 24754335-12d8-41f4-a44b-5f353e2a6ee2
av_lk = mean(ArDCA.loglikelihood(s, arnet) for s in aln)

# ╔═╡ fdce29fa-6f49-46f7-bc00-bd52c765eff1
av_E = mean(s -> DCATools.energy(potts, s), aln)

# ╔═╡ 0c440200-a26c-465d-813e-269f52d2a8e1
md"## Equilibrated init."

# ╔═╡ 633b097f-3676-48e3-9bf9-71777df880ce
AutoRegressiveASR.average_at_t

# ╔═╡ 52fdbaed-6081-41a1-a7b7-a56d975d5441
fam_picker

# ╔═╡ 0905421b-0dd7-4f31-aa95-319339f1b15b
eps()

# ╔═╡ de0eefad-5099-4835-96b2-69d1b285c749
md"# Functions"

# ╔═╡ e6eea28e-a8b4-49a1-96d2-f44dcc25315d
function my_fitexp(x, y)
	ymax = maximum(y) * (1.1)
	z = log.(1 .- y / ymax) # ~ exp(a*x)
	a = sum(x.*z)/sum(x.^2)
	return -1/a, ymax, t -> ymax * (1-exp(a*t))
end

# ╔═╡ 4423852f-b3c9-485c-9d93-2ae56933aa90
measures = AutoRegressiveASR.useful_callbacks()

# ╔═╡ df4c41c1-84ae-4c7c-9236-0de402318584
let
	function pubfig(fnt_size=30; kwargs...)
	    fnt_size = 30
	    PLOTS_DEFAULTS = Dict(
	        :markersize => 10,
	        :linewidth => 5,
	        :titlefontsize => fnt_size,
	        :guidefontsize => fnt_size,
	        :tickfontsize => fnt_size,
	        :legendfontsize => fnt_size,
	        :size => (1200,900),
	        :gridlinewidth => 0.5,
	        :gridalpha => 0.3,
	        :framestyle => :box,
	        :margin => 5mm,
	        # :bottom_margin => 5mm,
	        # :left_margin => 5mm,
	    )
	
	    for (k, x) in kwargs
	        PLOTS_DEFAULTS[k] = x
	    end
	
	    return PLOTS_DEFAULTS
	end
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ c73a1d58-1fbc-4912-ac92-73e03dc5569f
function logrange(x, y, l; add_zero=false, round_to_int=true) 
	R = @chain range(log(x), log(y); length=l) exp.(_)
	return collect(add_zero ? vcat(0, R) : R)
end

# ╔═╡ 7b5f3af1-1b0d-4c33-82cc-fd51437d3ebe
eq_measures = let
	tvals = logrange(.01, 10, 15; add_zero = true)
	nat_seqs = rand(aln, 100)
	vals = map(nat_seqs) do s0 
		AutoRegressiveASR.average_at_t(
			s0, measures, arnet, tvals, 50; ref_model=potts
		)
	end
end

# ╔═╡ 9cb086a9-4feb-4567-93ba-fff4b4b5a2d2
let p = plot()
	for T in eq_measures
		tvals = [x.t for x in T]
		E = [x.energy for x in T]
		plot!(tvals, E, label="", color = 1, line=(2, .5))
	end

	E = mean(T -> [x.energy for x in T], eq_measures)
	tvals = [x.t for x in first(eq_measures)]
	plot!(tvals, E, label="average", color=:black)
	hline!([av_E], line=(2, :black, :dash), label="")
	plot!(
		xscale=:log10, 
		xlim=(.01, maximum(tvals)),
		xlabel = "time",
		ylabel = "<E>",
		title = "Eq. init",
		xticks = [.01, .1, 1],
	)
end

# ╔═╡ 16996142-7f1e-4f90-8a01-0f17703f137e
let p = plot()
	for T in eq_measures
		tvals = [x.t for x in T]
		H = [x.hamming_to_init for x in T]
		plot!(tvals, H/L, label="", color = 1, line=(2, .5))
	end

	H = mean(T -> [x.hamming_to_init for x in T], eq_measures)/L
	tvals = [x.t for x in first(eq_measures)]
	plot!(tvals, H, label="average", color=:black)

	F = fitexp(tvals, H)
	plot!(tvals, F.(tvals), label="", line=(:black, :dash, 3))
	# @info F
	Teq = 1/F.b

	for (i, t) in enumerate([1.5, 2])
		h = F(t)
		plot!([t, t], [0, h]; line=(2, :dash), label="", color = i+1)
		plot!([1e-3, t], [h, h]; line = (2, :dash), label="", color = i+1)
	end
	h2 = F(2)
	
	plot!(
		xscale=:log10, 
		xlim=(.01, maximum(tvals)),
		xlabel = "time",
		ylabel = "Hamming to init",
		title = "Eq. init - Teq = $(Teq)",
		# xticks = [.1, 1, 10, 100, 1000, 10_000],
		legend = :topleft,
	)
end

# ╔═╡ 1e28b91f-ad87-45d1-8b4e-811240e7981f
let p = plot()
	H = mean(T -> [x.hamming_to_init for x in T], eq_measures)/L
	tvals = [x.t for x in first(eq_measures)]
	plot!(tvals, H, label="average", color=:black)

	F = fitexp(tvals, H)
	plot!(tvals, F.(tvals))
end

# ╔═╡ Cell order:
# ╠═ec42999a-edaf-11ee-21a4-8558cad7bc28
# ╠═5dcb047d-46b6-45d7-9780-a98c0055eb7a
# ╠═7df58d8b-95a4-4a2c-ad5f-a1747797598c
# ╠═dbdf3b60-8579-472d-8236-7005788d3ceb
# ╠═559526a0-8933-4508-81e4-6a61d6e1d229
# ╠═eadcdce5-dd16-409e-8e7b-d70edc2eb380
# ╠═57d52ec7-6085-4c44-aeef-0bd48baa156a
# ╠═36e7cb69-22ac-4c09-bdf3-37c217b03d29
# ╠═92b6f8a2-b01a-4a9e-9167-4a4a80b14966
# ╠═24754335-12d8-41f4-a44b-5f353e2a6ee2
# ╠═fdce29fa-6f49-46f7-bc00-bd52c765eff1
# ╠═0c440200-a26c-465d-813e-269f52d2a8e1
# ╠═633b097f-3676-48e3-9bf9-71777df880ce
# ╠═7b5f3af1-1b0d-4c33-82cc-fd51437d3ebe
# ╠═9cb086a9-4feb-4567-93ba-fff4b4b5a2d2
# ╠═52fdbaed-6081-41a1-a7b7-a56d975d5441
# ╠═16996142-7f1e-4f90-8a01-0f17703f137e
# ╠═0905421b-0dd7-4f31-aa95-319339f1b15b
# ╠═1e28b91f-ad87-45d1-8b4e-811240e7981f
# ╟─de0eefad-5099-4835-96b2-69d1b285c749
# ╠═e6eea28e-a8b4-49a1-96d2-f44dcc25315d
# ╠═4423852f-b3c9-485c-9d93-2ae56933aa90
# ╠═df4c41c1-84ae-4c7c-9236-0de402318584
# ╠═c73a1d58-1fbc-4912-ac92-73e03dc5569f
