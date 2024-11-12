### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ 25ed7ab2-085d-11ef-3f81-25b122d92163
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

# ╔═╡ 35623ac1-4a5a-441d-a37d-5f8f41c9159d
begin
	local plot_defaults = joinpath(homedir(), ".julia/config/plot_defaults.jl")
	if isfile(plot_defaults)
		include(plot_defaults)
	end
end

# ╔═╡ a29c5186-57c3-4713-a17c-3a5b79ab3dd9
if isdefined(@__MODULE__, :pubfig)
	plt_defaults = pubfig(20)
	Plots.default(; plt_defaults...)
end

# ╔═╡ c1328117-d0a1-4482-8cc9-1e5ebd621f05
# pubfig(20)

# ╔═╡ 4afe6959-60ee-4463-ad56-56e022588ac0
md"## Folders"

# ╔═╡ 11db0774-812b-4107-98c4-cb3c17e2c365
folder_list = let
    F = vcat(
    	# readdir(datadir("simulated/arnet_yule"); join=true),
    	readdir(datadir("simulated/potts_yule"); join=true),
    )
    filter(f -> occursin(r"PF\d\d\d\d\d", basename(f)), F)
end

# ╔═╡ 297331e9-c210-447b-9d55-ba20880170c9
folder_picker = @bind folder_full Select(folder_list)

# ╔═╡ 22b784c2-1e0d-43ed-a7b3-f4de62ba06d6
folder = basename(folder_full)

# ╔═╡ 8d8207d3-2ea6-482a-a646-bc1e6bcb6f74
fam_main = "PF00072"

# ╔═╡ 83321a10-4baf-43a7-be13-d67d86873ac4
md"# Individual plots"

# ╔═╡ 00b47571-ce7d-493f-a40f-37c323fe0157
includet(
    projectdir("notebooks/article_figures/base_notebooks/hamming_distances_as_func.jl")
)
# HAM = pluto_ingredients(projectdir(
# 	"notebooks/article_figures/base_notebooks/hamming_distances_as_func.jl"
# ))

# ╔═╡ a0a5e374-194c-4872-865b-abe367036859
iqtree_strategies = ["base", "C10", "C60"]

# ╔═╡ a8a77324-8462-4afe-bfe3-01122bf0048a

# ╔═╡ aa3ffef5-4a52-4e3e-a0f7-beb94939485e
function make_plot(folder, global_title; panel=true, kwargs...)
    ham_plots = hamming_distance_plots(folder; kwargs...)
	# ham_plots = HAM.hamming_distance_plots(folder)
	p1, p2, p3 = ham_plots

	plot!(p1, title = "MAP reconstruction")
	plot!(p2, title = "MAP reconstruction - no gaps")
	plot!(p3, title = "Posterior sampling - no gaps")

	split(folder, "_")[1]
    return if panel
    	plot(
    		p1, p2, p3;
    		layout = grid(1, 3),
            size=(2700,900),
    		left_margin=15mm,
            bottom_margin=15mm,
    		dpi = 300,
    		plot_title = global_title,
            plot_titlefontsize = 26,
    	)
    else
        p1, p2, p3
	end
end

# ╔═╡ afd4850a-1ecb-4952-8268-fc143a59bbc9
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
	    savename = "hamming_to_real_$(evolver)_$(fam)_iqtree-$(iqtree_strategy)"
	    p = make_plot(folder, title; iqtree_strategy)
	    # savefig(p, projectdir(savedir, savename*".png"))
        plot!(gridalpha=0.05)
        savefig(p, projectdir(savedir, savename*".pdf"))
    end
end

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
plts = Dict()
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


    for iqtree_strategy in iqtree_strategies
        @info iqtree_strategy
        savedir = projectdir("notes/article/figures/SI/")
        savename = "hamming_to_real_$(evolver)_extrafams_iqtree-$(iqtree_strategy).pdf"
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
            ps = make_plot(folder, fam; iqtree_strategy, panel=false)
            foreach(p -> plot!(p, bottom_margin=15mm, left_margin=10mm, top_margin=-2mm), ps)
            title = plot(
                title = fam,
                grid = false,
                showaxis = false,
                bottom_margin = -15mm,
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


# ╔═╡ Cell order:
# ╠═25ed7ab2-085d-11ef-3f81-25b122d92163
# ╠═35623ac1-4a5a-441d-a37d-5f8f41c9159d
# ╠═a29c5186-57c3-4713-a17c-3a5b79ab3dd9
# ╠═c1328117-d0a1-4482-8cc9-1e5ebd621f05
# ╠═4afe6959-60ee-4463-ad56-56e022588ac0
# ╠═11db0774-812b-4107-98c4-cb3c17e2c365
# ╠═22b784c2-1e0d-43ed-a7b3-f4de62ba06d6
# ╠═297331e9-c210-447b-9d55-ba20880170c9
# ╠═8d8207d3-2ea6-482a-a646-bc1e6bcb6f74
# ╠═83321a10-4baf-43a7-be13-d67d86873ac4
# ╠═00b47571-ce7d-493f-a40f-37c323fe0157
# ╠═a0a5e374-194c-4872-865b-abe367036859
# ╠═a8a77324-8462-4afe-bfe3-01122bf0048a
# ╠═afd4850a-1ecb-4952-8268-fc143a59bbc9
# ╠═aa3ffef5-4a52-4e3e-a0f7-beb94939485e
# ╠═f85148ca-d1f0-481c-b182-7f1bbf483071
# ╠═b674036b-be81-4f90-8635-7f5e12d9a7c2
# ╠═c871f90f-667f-4b7c-a75d-05adb8c2d7ad
# ╠═20d758be-c468-4e0a-8ee1-c6e47eaa2af0
