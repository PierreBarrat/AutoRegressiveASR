### A Pluto.jl notebook ###
# v0.19.42

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

# ╔═╡ 3d127b7a-c983-11ee-2228-93a82a0887ef
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")
	using AutoRegressiveASR
	using Chain
	using DCATools
	using EasyFit
	using Measures
	using Plots
	using PlutoUI
	using StatsBase
end

# ╔═╡ 859d9856-3f9d-4ed5-8854-08160aa1d63d
families = let
	FAM = pluto_ingredients(scriptsdir("families.jl"))
	FAM.families
end

# ╔═╡ 6eb8f929-49b5-4bcb-ace7-14684c1114cb
fam_picker = @bind family Select(map(reverse, collect(pairs(families))))

# ╔═╡ 33c40657-310d-4e32-8c47-6ef4b7ebe037
family["potts"]

# ╔═╡ ae0da5bb-6a48-4112-af93-bbe3f807cef0
potts = begin
	potts_file = family["potts"]
	DCAGraph(potts_file)
end

# ╔═╡ ce9fdcdc-9b6d-42c6-9a3c-142841cd05b0
aln = read_msa(family["sample_potts_eq"]);

# ╔═╡ 393decf5-4626-46e4-8daf-913cdc2ad690
begin
	q = potts.q
	L = potts.L
end

# ╔═╡ 84581806-cbed-4d27-837f-5b71aefb5b00
pw_distance_training, av_pw_distance_training = begin
	X = DCATools.pw_hamming_distance(aln; step=100)
	ecdf(X), mean(X)
end

# ╔═╡ b3d252a0-3b96-41f2-957e-e49c4d4cbffe
md"## Random init."

# ╔═╡ 85e96372-aac2-4de0-a92f-a64bf720f7de
fam_picker

# ╔═╡ 266cfb49-d735-47c2-873d-135226526414
md"## Equilibrated init."

# ╔═╡ 3e396945-58c1-4d8b-bea4-e6cacd216ca1


# ╔═╡ 19c4ca93-8068-4f89-86e8-87f084d5b127


# ╔═╡ 61e239ef-3c8d-4788-bf5c-0f9b9ed3775c
function my_fitexp(x, y)
	ymax = maximum(y) * (1.1)
	z = log.(1 .- y / ymax) # ~ exp(a*x)
	a = sum(x.*z)/sum(x.^2)
	return -1/a, ymax, t -> ymax * (1-exp(a*t))
end

# ╔═╡ 3dc16393-5ec9-4b05-a785-666100e13dc6
md"## Hamming distance distribution"

# ╔═╡ 2aaa227a-f643-45ef-950e-fc3006951add
md"# Functions"

# ╔═╡ 043e9a03-d25a-400c-a4a6-4e28bf0648ae
measures = AutoRegressiveASR.useful_callbacks()

# ╔═╡ 711776ad-0cd0-4119-84f3-defae8ab740a
measures

# ╔═╡ 6211b906-f089-4e7e-9c86-bda6ca9870f9
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

# ╔═╡ e33b6fcb-3a4a-478a-85f8-7b8f9caf5aad
function logrange(x, y, l; add_zero=false, round_to_int=true) 
	R = @chain range(log(x), log(y); length=l) exp.(_) round.(Int, _)
	return collect(add_zero ? vcat(0, R) : R)
end

# ╔═╡ 6e56e613-1658-410b-a125-308b386e09e5
RI_measures = let
	tvals = logrange(1, 100_000, 12; add_zero = true)
	random_inits = rand(1:q, L, 5)
	vals = map(eachcol(random_inits)) do s0 
		AutoRegressiveASR.average_at_t(s0, measures, potts, tvals, 25)
	end
end

# ╔═╡ be0fd821-a066-4d48-9a17-30533a38efed
let p = plot()
	for T in RI_measures
		tvals = [x.t for x in T]
		E = [x.energy for x in T]
		plot!(tvals, E, label="", color = 1, line=(2, .5))
	end

	E = mean(T -> [x.energy for x in T], RI_measures)
	tvals = [x.t for x in first(RI_measures)]
	plot!(tvals, E, label="average", color=:black)
	
	plot!(
		xscale=:log10, 
		xlim=(1, 2*maximum(tvals)),
		xlabel = "mcmc steps",
		ylabel = "<E>",
		title = "Random init",
		xticks = [1, 100, 10_000, ]
	)
end

# ╔═╡ 809fabf1-5807-4f09-b02e-85bd5654e874
let p = plot()
	for T in RI_measures
		tvals = [x.t for x in T]
		H = [x.hamming_to_init for x in T]
		plot!(tvals/L, H/L, label="", color = 1, line=(2, .5))
	end

	H = mean(T -> [x.hamming_to_init for x in T], RI_measures)
	tvals = [x.t for x in first(RI_measures)]
	plot!(tvals/L, H/L, label="average", color=:black)

	hinf = maximum(H)/L
	hline!([0.75 * hinf], label = "", line = (:black, :dash, 1))
	tt = +(
		tvals[findlast(<(0.75*hinf), H/L)]/L,
		tvals[findfirst(>(0.75*hinf), H/L)]/L,
	)/2
	plot!([tt, tt], [0, 0.75*hinf], label="")
	@info tt
	
	plot!(
		xscale=:log10, 
		xlim=(1/L, 2*maximum(tvals)/L),
		xlabel = "mcmc sweeps",
		ylabel = "Hamming to init",
		title = "Random init",
		xticks = [.1, 1, 10, 100, 10_000, ],
		legend = :bottomright,
	)
end

# ╔═╡ e55261cf-108d-4a36-beb8-cdb0e16fbafd
let p = plot()
	for T in RI_measures
		tvals = [x.t for x in T]/L
		H = [x.pw_hamming for x in T]
		plot!(tvals, H/L, label="", color = 1, line=(2, .5))
	end

	H = mean(T -> [x.pw_hamming for x in T], RI_measures)
	tvals = [x.t for x in first(RI_measures)]/L
	plot!(tvals, H/L, label="average", color=:black)
	
	plot!(
		xscale=:log10, 
		xlim=(1/L, 2*maximum(tvals)),
		xlabel = "mcmc sweeps",
		ylabel = "Self-Hamming",
		title = "Random init",
		xticks = [.1, 1, 10, 100, 10_000, ],
		legend = :bottomright,
	)
end

# ╔═╡ 01765a1f-c8ba-454b-8df9-7277a115d29a
eq_measures = let
	tvals = logrange(1, 50_000, 12; add_zero = true)
	nat_seqs = rand(aln, 25)
	vals = map(nat_seqs) do s0 
		AutoRegressiveASR.average_at_t(s0, measures, potts, tvals, 50)
	end
end

# ╔═╡ 8c5ebcde-2cc0-42af-8e48-b02daf15b9f9
let p = plot()
	for T in eq_measures
		tvals = [x.t for x in T]
		E = [x.energy for x in T]
		plot!(tvals/L, E, label="", color = 1, line=(2, .5))
	end

	E = mean(T -> [x.energy for x in T], eq_measures)
	tvals = [x.t for x in first(eq_measures)]
	plot!(tvals/L, E, label="average", color=:black)
	
	plot!(
		xscale=:log10, 
		xlim=(1, 2*maximum(tvals/L)),
		xlabel = "mcmc sweeps",
		ylabel = "<E>",
		title = "Eq. init",
		xticks = [1, 10, 100, 1_000],
	)
end

# ╔═╡ 06ee81ad-ce1d-4254-bd8f-c6df2bce1176
let p = plot()
	for T in eq_measures
		tvals = [x.t for x in T]
		H = [x.hamming_to_init for x in T]
		plot!(tvals/L, H/L, label="", color = 1, line=(2, .5))
	end

	H = mean(T -> [x.hamming_to_init for x in T], eq_measures)
	tvals = [x.t for x in first(eq_measures)]
	plot!(tvals/L, H/L, label="average", color=:black)

	# Teq, hmax, expfit = my_fitexp(tvals/L, H/L)
	# @info Teq, hmax
	# plot!(tvals/L, expfit.(tvals/L), line = (:black, :dash), label="fit")
	# Teq = round(expfit.b, sigdigits=2)
	
	plot!(
		xscale=:log10, 
		xlim=(1/L, 2*maximum(tvals)/L),
		xlabel = "mcmc sweeps",
		ylabel = "Hamming to init",
		# title = "Eq. init - Teq = $(Teq) sweeps",
		xticks = [.1, 1, 10, 100, 1000, 10_000],
		legend = :bottomright,
	)
end

# ╔═╡ 33842c68-205e-41ac-8368-57258480e236
let p = plot()
	for T in eq_measures
		tvals = [x.t for x in T]/L
		H = [x.pw_hamming for x in T]
		plot!(tvals, H/L, label="", color = 1, line=(2, .5))
	end

	H = mean(T -> [x.pw_hamming for x in T], eq_measures)
	tvals = [x.t for x in first(eq_measures)]/L
	plot!(tvals, H/L, label="average", color=:black)
	
	plot!(
		xscale=:log10, 
		xlim=(1/L, 2*maximum(tvals)),
		xlabel = "mcmc sweeps",
		ylabel = "Self-Hamming",
		title = "Eq. init",
		xticks = [.1, 1, 10, 100, 10_000, ],
		legend = :bottomright,
	)
end

# ╔═╡ bac4fb4a-5409-43dc-a8e2-78b3ffef721d
hamming_dist_v_t = let
	s0 = rand(aln, 1)[1]
	tvals=  logrange(1, 320, 6) * L
	samples = AutoRegressiveASR.propagate(s0, tvals, potts, 100)
	map(zip(tvals, samples)) do (t, S)
		(t=t, H_cdf = ecdf(DCATools.pw_hamming_distance(S)))
	end
end

# ╔═╡ b4161cde-73df-459d-b892-3c32a47885c9
let p = plot()
	pal = palette(:roma, length(hamming_dist_v_t))
	xvals = 0:.01:1
	plot!(
		xvals, pw_distance_training.(xvals);
		line = (:black, :dash), label="natural",
	)
	for (i, (t, d)) in enumerate(hamming_dist_v_t)
		plot!(
			xvals, d.(xvals);
			label="t=$(round(t/L; sigdigits=2))", color = pal[i], line=(4)
		)
	end
	plot!(
		xlabel = "pw hamming distance",
		ylabel = "",
		title = "cdf of pw hamming distance \n time in mcmc sweeps",
		legend = :bottomleft,
		size = (1200,1200),
	)
end

# ╔═╡ Cell order:
# ╠═3d127b7a-c983-11ee-2228-93a82a0887ef
# ╠═859d9856-3f9d-4ed5-8854-08160aa1d63d
# ╠═6eb8f929-49b5-4bcb-ace7-14684c1114cb
# ╠═33c40657-310d-4e32-8c47-6ef4b7ebe037
# ╠═ae0da5bb-6a48-4112-af93-bbe3f807cef0
# ╠═ce9fdcdc-9b6d-42c6-9a3c-142841cd05b0
# ╠═393decf5-4626-46e4-8daf-913cdc2ad690
# ╠═84581806-cbed-4d27-837f-5b71aefb5b00
# ╟─b3d252a0-3b96-41f2-957e-e49c4d4cbffe
# ╠═6e56e613-1658-410b-a125-308b386e09e5
# ╟─be0fd821-a066-4d48-9a17-30533a38efed
# ╠═85e96372-aac2-4de0-a92f-a64bf720f7de
# ╠═809fabf1-5807-4f09-b02e-85bd5654e874
# ╟─e55261cf-108d-4a36-beb8-cdb0e16fbafd
# ╟─266cfb49-d735-47c2-873d-135226526414
# ╠═01765a1f-c8ba-454b-8df9-7277a115d29a
# ╠═711776ad-0cd0-4119-84f3-defae8ab740a
# ╠═3e396945-58c1-4d8b-bea4-e6cacd216ca1
# ╟─8c5ebcde-2cc0-42af-8e48-b02daf15b9f9
# ╠═06ee81ad-ce1d-4254-bd8f-c6df2bce1176
# ╠═19c4ca93-8068-4f89-86e8-87f084d5b127
# ╠═61e239ef-3c8d-4788-bf5c-0f9b9ed3775c
# ╟─33842c68-205e-41ac-8368-57258480e236
# ╟─3dc16393-5ec9-4b05-a785-666100e13dc6
# ╠═bac4fb4a-5409-43dc-a8e2-78b3ffef721d
# ╟─b4161cde-73df-459d-b892-3c32a47885c9
# ╟─2aaa227a-f643-45ef-950e-fc3006951add
# ╠═043e9a03-d25a-400c-a4a6-4e28bf0648ae
# ╠═6211b906-f089-4e7e-9c86-bda6ca9870f9
# ╠═e33b6fcb-3a4a-478a-85f8-7b8f9caf5aad
