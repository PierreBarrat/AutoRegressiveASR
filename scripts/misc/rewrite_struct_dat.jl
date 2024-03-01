### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 2d693f40-d722-11ee-1095-77ddc1be3bd5
begin
	using DrWatson
	quickactivate(@__DIR__, "AutoRegressiveASR")
	using CSV
	using DataFrames
	using DataFramesMeta
end

# ╔═╡ 72df1702-9c47-4a59-a5f3-72c852995e48
fam = "PF00595"
struct_file = datadir("Models/$fam/$(fam)_struct.dat")

function rewrite_struct_file(struct_file)
    dat_init = CSV.File(struct_file; header=false,) |> DataFrame
    dat_out = @orderby DataFrame(
    	:i => Int.(dat_init[:, 1]),
    	:j => Int.(dat_init[:, 2]),
    	:distance_A => dat_init[:,4],
    ) :i :j
    outfile = joinpath(dirname(struct_file), splitext(struct_file)[1] * ".csv")
    CSV.write(outfile, dat_out, delim=' ')
end



