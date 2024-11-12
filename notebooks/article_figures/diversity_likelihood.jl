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

# ╔═╡ 56492eb5-a7aa-42c1-ade3-90cd3cccdbf2
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")
end
begin
	using AutoRegressiveASR
	using Chain
	using Measures
	using Plots
	using PlutoUI
end

# ╔═╡ 7e7d1989-77dd-49ad-86a1-fdcccf5ed8f4
begin
	local plot_defaults = joinpath(homedir(), ".julia/config/plot_defaults.jl")
	if isfile(plot_defaults)
		include(plot_defaults)
	end
end

# ╔═╡ 026c32c2-86de-42b8-a6a0-865b155094bc
if isdefined(@__MODULE__, :pubfig)
	plt_defaults = pubfig(22)
	Plots.default(; plt_defaults...)
end

# ╔═╡ d3dddbb2-f1cf-43bd-b936-aa9281ff4605


# ╔═╡ 13fdbc6d-4497-4394-bd9b-ec662f3f625e
md"## Folders"

# ╔═╡ 4a9db40c-f0ce-49ac-b642-5cb432ec9f01
folder_list = let
    F = vcat(
        readdir(datadir("simulated/arnet_yule"); join=true),
        # readdir(datadir("simulated/potts_yule"); join=true),
    )
    filter(f -> occursin(r"PF\d\d\d\d\d", basename(f)), F)
end

# ╔═╡ 3c20b9aa-f4cd-4429-ac91-4152cd22a9ce

# ╔═╡ 568efcdc-7cb3-4229-96d3-df2d08b33a40

# ╔═╡ 307ec663-ac37-4571-b55d-0ce223143286
fam_main = "PF00072"

# ╔═╡ e7f68da4-6a02-4ef9-a392-02ccbd5e4912
md"# Figures"

# ╔═╡ 6a91f6c3-0e3d-456d-95df-1b18ef5b3e26
includet(projectdir("notebooks/article_figures/base_notebooks/diversity_as_func.jl"))

# ╔═╡ a2ef67da-a479-4751-990e-dd140cbc4d8c
includet(
    projectdir(
	    "notebooks/article_figures/base_notebooks/likelihood_and_consensus_as_func.jl"
    )
)

# ╔═╡ 057c0914-adbc-4dad-a16d-8697d8a9cef4
iqtree_strategies = ["base", "C10", "C60"]

# ╔═╡ b041d20c-bf5c-4563-b473-87fede4a7b94
function make_plot(folder, global_title; panel=true, kwargs...)
	self_hamming, _ = diversity_plots(folder; kwargs...)
	lk, cons = likelihood_and_hamming_consensus(folder; kwargs...)

	# plot!(lk, left_margin=20mm)
	# plot!(self_hamming, left_margin=15mm)
	# plot!(cons, left_margin=20mm)

	return if panel
        plot(
    		self_hamming, cons, lk;
    		layout = (@layout [a{0.31w} b{0.31w} c{0.38w}]),
    		size = (2400, 800),
    		bottom_margin = 15mm, left_margin = 25mm,
    		dpi = 300,
    		plot_title = global_title,
            plot_titlefontsize = 26,
    	)
    else
        self_hamming, cons, lk
    end
end

# ╔═╡ 5102c585-4998-4d79-b649-2cdafafedb38
for folder in folder_list
	fam = split(basename(folder), "_")[1]
    fam != fam_main && continue
	evolver = @chain folder begin
		splitpath
		_[end-1]
		split("_")
		getindex(1)
	end

	savedir = projectdir(
		"notes/article/figures/",
		(fam == fam_main) ? "." : "SI"
	)
	title = fam == fam_main ? "" : fam

    for iqtree_strategy in iqtree_strategies
    	savename = "diversity_likelihood_$(evolver)_$(fam)_iqtree-$(iqtree_strategy)"
    	p = make_plot(folder, title; iqtree_strategy)
    	# evolver == "arnet" && plot!(p, xlim = (-0.025, 2.025))
    	# savefig(p, projectdir(savedir, savename*".png"))
        plot!(gridalpha=0.05)
        savefig(p, projectdir(savedir, savename*".pdf"))
    end
end

#===============================================================#
###################### Extra family panels ######################
#===============================================================#

sim_parameters, families = let
    params_fams = map(folder_list) do f
       fam = match(r"PF\d\d\d\d\d", f).match |> string
       (split(f, fam), fam)
   end
   params = unique([x[1] for x in params_fams])
   families = unique([x[2] for x in params_fams])

   params, families
end

plts = Dict()
if isdefined(@__MODULE__, :pubfig)
    plt_defaults = pubfig(14)
    Plots.default(; plt_defaults...)
end
for prms in sim_parameters
    # prms[1] * fam * prms[2] is an element of folder_list (abs path)
    # prms[1] is the root path, with arnet/potts in it
    # prms[2] has the rest of the info
    # /home/.../arnet_yule/PFXXX_nleaves=...
    @info prms
    if !all(fam -> isdir(prms[1] * fam * prms[2]), families)
        continue
    end
    evolver = @chain splitpath(prms[1])[end] split("_") getindex(1)
    evolver == "potts" && continue


    for iqtree_strategy in iqtree_strategies
        @info iqtree_strategy
        savedir = projectdir("notes/article/figures/SI/")
        savename = "diversity_likelihood_$(evolver)_extrafams_iqtree-$(iqtree_strategy).pdf"
        savepath = joinpath(savedir, savename)
        if isfile(savepath)
            @error "$savepath already exists - skipping"
            continue
        end
        @info "Save in $(savepath)"
        # fam specific panel
        fam_plts = []
        for fam in filter(!=(fam_main), families)
            folder = prms[1] * fam * prms[2]
            ps = make_plot(folder, ""; iqtree_strategy, panel=false)
            foreach(p -> plot!(p, bottom_margin=15mm, left_margin=10mm), ps)
            title = plot(
                title = fam,
                grid = false,
                showaxis = false,
                bottom_margin = -23mm,
                titlefontsize=24,
            )
            p = plot(
                title, ps...;
                layout=@layout([A{0.01h}; [B C D]]),
                size = (1800, 900),
                dpi=300,
            )
            push!(fam_plts, p)
        end
        panel = plot(fam_plts...; layout=grid(3,1), size=(1600, 1600), gridalpha=0.05)
        savefig(panel, savepath)
        plts[prms, iqtree_strategy] = fam_plts
    end
end


# ╔═╡ fda52b93-7e1c-4e3d-a346-181bbc7ff815
# p = make_plot(folder_full, split(basename(folder_full), "_")[1])

# ╔═╡ 0c3c40b5-fd1a-4034-b5f1-88bab4db081e

# ╔═╡ 07436b06-dd3e-48ad-92fe-b790ff9c66ab


# ╔═╡ Cell order:
# ╠═56492eb5-a7aa-42c1-ade3-90cd3cccdbf2
# ╠═7e7d1989-77dd-49ad-86a1-fdcccf5ed8f4
# ╠═026c32c2-86de-42b8-a6a0-865b155094bc
# ╠═d3dddbb2-f1cf-43bd-b936-aa9281ff4605
# ╠═13fdbc6d-4497-4394-bd9b-ec662f3f625e
# ╠═4a9db40c-f0ce-49ac-b642-5cb432ec9f01
# ╠═568efcdc-7cb3-4229-96d3-df2d08b33a40
# ╠═3c20b9aa-f4cd-4429-ac91-4152cd22a9ce
# ╠═307ec663-ac37-4571-b55d-0ce223143286
# ╠═e7f68da4-6a02-4ef9-a392-02ccbd5e4912
# ╠═6a91f6c3-0e3d-456d-95df-1b18ef5b3e26
# ╠═057c0914-adbc-4dad-a16d-8697d8a9cef4
# ╠═a2ef67da-a479-4751-990e-dd140cbc4d8c
# ╠═5102c585-4998-4d79-b649-2cdafafedb38
# ╠═b041d20c-bf5c-4563-b473-87fede4a7b94
# ╠═fda52b93-7e1c-4e3d-a346-181bbc7ff815
# ╠═0c3c40b5-fd1a-4034-b5f1-88bab4db081e
# ╠═07436b06-dd3e-48ad-92fe-b790ff9c66ab
