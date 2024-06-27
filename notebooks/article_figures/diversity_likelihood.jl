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
pubfig(20)

# ╔═╡ 13fdbc6d-4497-4394-bd9b-ec662f3f625e
md"## Folders"

# ╔═╡ 4a9db40c-f0ce-49ac-b642-5cb432ec9f01
folder_list = vcat(
	readdir(datadir("simulated/arnet_yule"); join=true),
	readdir(datadir("simulated/potts_yule"); join=true),
)

# ╔═╡ 3c20b9aa-f4cd-4429-ac91-4152cd22a9ce
folder_picker = @bind folder_full Select(folder_list)

# ╔═╡ 568efcdc-7cb3-4229-96d3-df2d08b33a40
folder = basename(folder_full)

# ╔═╡ 307ec663-ac37-4571-b55d-0ce223143286
fam_main = "PF00072"

# ╔═╡ e7f68da4-6a02-4ef9-a392-02ccbd5e4912
md"# Figures"

# ╔═╡ 6a91f6c3-0e3d-456d-95df-1b18ef5b3e26
DIV = pluto_ingredients(projectdir(
	"notebooks/article_figures/base_notebooks/diversity_as_func.jl"
))

# ╔═╡ 057c0914-adbc-4dad-a16d-8697d8a9cef4
folder_full

# ╔═╡ a2ef67da-a479-4751-990e-dd140cbc4d8c
LK = pluto_ingredients(projectdir(
	"notebooks/article_figures/base_notebooks/likelihood_and_consensus_as_func.jl"
))

# ╔═╡ b041d20c-bf5c-4563-b473-87fede4a7b94
function make_plot(folder, global_title)
	self_hamming, _ = DIV.diversity_plots(folder)
	lk, cons = LK.likelihood_and_hamming_consensus(folder)

	# plot!(lk, left_margin=20mm)
	# plot!(self_hamming, left_margin=15mm)
	# plot!(cons, left_margin=20mm)

	p = plot(
		self_hamming, cons, lk;
		layout = (@layout [a{0.31w} b{0.31w} c{0.38w}]),
		size = (2400, 800), 
		bottom_margin = 15mm, left_margin = 25mm,
		dpi = 300,
		plot_title = global_title,
	)
	plot!(plot_titlefontsize = 26)
	
	p
end

# ╔═╡ 5102c585-4998-4d79-b649-2cdafafedb38
for folder in folder_list
	fam = split(basename(folder), "_")[1]
	method = @chain folder begin
		splitpath
		_[end-1]
		split("_")
		getindex(1)
	end

	savename = "diversity_likelihood_$(method)_$(fam).png"
	savedir = projectdir(
		"notes/article/figures/",
		(fam == fam_main) ? "." : "SI"
	)
	title = fam == fam_main ? "" : fam

	p = make_plot(folder, title)
	# method == "arnet" && plot!(p, xlim = (-0.025, 2.025))
	savefig(p, projectdir(savedir, savename))
	p
end

# ╔═╡ fda52b93-7e1c-4e3d-a346-181bbc7ff815
p = make_plot(folder_full, split(basename(folder_full), "_")[1])

# ╔═╡ 0c3c40b5-fd1a-4034-b5f1-88bab4db081e
savefig(p, "/home/pierrebc/Bureau/plot.png")

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
