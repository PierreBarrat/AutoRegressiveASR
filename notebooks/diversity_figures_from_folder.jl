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

# ╔═╡ 988a4668-d718-11ee-1816-b73ece48c950
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")

	using AutoRegressiveASR
	using DataFrames
	using DataFramesMeta
	using PlutoUI
	using StatsBase
	using StatsPlots
end

# ╔═╡ 3a06c22c-bbfe-41c1-9b68-26482fc1cc34
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 79db37a3-ed65-41f0-9e95-275bf7fac3a2
begin
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 7990e53e-27cc-4d31-9d1b-c2c172d20d5c
_fs = @bind folder_full Select(readdir(datadir("simulated/potts_yule"); join=true))

# ╔═╡ 7c3a2726-4925-4ed6-b8fc-dcc02f8f466d
path_folder, data_folder, prefix = let
	p, d = dirname(folder_full), basename(folder_full)
	prefix = split(d, "_")[1]
	p, d, prefix
end

# ╔═╡ 1d083c52-ccba-44fc-99ea-035ecefc7abf
prefix

# ╔═╡ 70a448e6-7e0f-4b4e-be30-05a5d4b30796
md"# Functions"

# ╔═╡ 8a4cab3c-892d-4d25-b72a-55f0f2da9be8
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

# ╔═╡ 13fef31f-abaa-4e94-917b-943b1b92901d
DIV = ingredients(scriptsdir("figures_and_results/diversity_functions.jl"))

# ╔═╡ 077bd775-60ce-41b3-bdaa-771b7784788d
data, filename = produce_or_load(
   Dict("basefolder" => folder_full, "out" => "diversity_data.jld2"); 
   filename = x -> joinpath(x["basefolder"], x["out"]),
   suffix = "",
) do config
   data = DIV.diversity_data(config["basefolder"], config["out"])
end;

# ╔═╡ 66c8ab09-20f0-4a22-83ec-291e424bb518
let p = plot()
	@df data["iqtree"] scatter!(:depth, :av_self_hamming)
	@df data["ardca"] scatter!(:depth, :av_self_hamming)
end

# ╔═╡ Cell order:
# ╠═988a4668-d718-11ee-1816-b73ece48c950
# ╠═3a06c22c-bbfe-41c1-9b68-26482fc1cc34
# ╠═79db37a3-ed65-41f0-9e95-275bf7fac3a2
# ╠═13fef31f-abaa-4e94-917b-943b1b92901d
# ╠═7990e53e-27cc-4d31-9d1b-c2c172d20d5c
# ╠═7c3a2726-4925-4ed6-b8fc-dcc02f8f466d
# ╠═1d083c52-ccba-44fc-99ea-035ecefc7abf
# ╠═077bd775-60ce-41b3-bdaa-771b7784788d
# ╠═66c8ab09-20f0-4a22-83ec-291e424bb518
# ╠═70a448e6-7e0f-4b4e-be30-05a5d4b30796
# ╠═8a4cab3c-892d-4d25-b72a-55f0f2da9be8
