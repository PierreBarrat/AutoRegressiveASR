using DrWatson
@quickactivate "AutoRegressiveASR"

if isempty(ARGS)
    println("""
        Usage: `infer_arnets.jl sequences.fasta outfolder`
        Results will be written in `dir/outfolder`, \
        with `dir` being the directory of `sequences.fasta`.
    """)
    exit()
end

using AutoRegressiveASR
using ArDCA
using Chain
using DCATools # necessary for "-A...Y" mapping
using DelimitedFiles
using JLD2
using JSON3

## Dealing with folders etc...
fastafile = ARGS[1] |> abspath
base_folder = dirname(fastafile) |> abspath
outfolder = joinpath(base_folder, ARGS[2]) |> abspath

@info "Learning ArDCA on sequences $fastafile. Results written in $outfolder."

weight_file = @chain begin
    splitext(fastafile)
    _[1] * "_weights" * _[2]
    abspath
end
@info "Looking for weights in $(weight_file)"
weight_file = isfile(weight_file) ? weight_file : nothing
isnothing(weight_file) && @info "No weights found, will compute them."

## Parameters for learning log
timestamp = now_string(; minute=true)
script = relpath(abspath(@__FILE__), projectdir())
parameters = @dict(
    fastafile, weight_file, timestamp, script
)

## Reading data and learning
S = read_msa(fastafile) |> unique
w = isnothing(weight_file) ? computeweights(S) : vec(readdlm(weight_file))
w = w/sum(w)

λvals = [1e-2, 1e-3, 1e-4]

for λ in λvals
    @info "Infering arnet for $fastafile (weights $weight_file) and regularization $λ"

    lambdaJ = λ
    lambdaH = λ/10

    dest = joinpath(outfolder, "arnet_lJ$(lambdaJ)_lH$(lambdaH).jld2")
    @info "Output will be saved at $dest"

    arnet, arvar = ardca(S.dat, w; lambdaJ, lambdaH)
    tagsave(dest, Dict("arnet" => arnet, "arvar" => arvar))
end

open(io -> JSON3.pretty(io, parameters), joinpath(outfolder, "parameters.json"), "w")
