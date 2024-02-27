using DrWatson
@quickactivate "AutoRegressiveASR"

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--outfolder"
        help = "Defaults to data/simulated/potts_yule/somename"
        arg_type = String
        default = datadir("simulated", "potts_yule")
    "--sample_potts_eq"
        help = "fasta file with a sample from the teacher potts model."
        arg_type = String
        default = ""
    "--nsweeps"
        help = "desired number of MCMC sweeps between root and leaves"
        arg_type = Int
        default = 1
    "--prefix"
        help = "prefix for output folder"
        arg_type = String
        default = ""
    "--add_outgroup"
        help = "Add outgroup for rooting"
        action = :store_true
    "--ntrees"
        help = "Number of trees to simulate"
        arg_type = Int
        default = 1
    "potts"
        help = "file containing the potts model"
        arg_type = String
        required = true
end
parsed_args = parse_args(ARGS, s)

using BackwardCoalescent
using DCATools
using JSON3

## Parameters
# Potts model
potts_file = parsed_args["potts"]
sample_potts_file = parsed_args["sample_potts_eq"]

# Genealogy
nleaves = 50 # number of leaves
nsweeps = parsed_args["nsweeps"] # desired number of MCMC sweeps between root and leaves
L = potts.L
b = log(n)/nsweeps/potts.L
outgroup = parsed_args["add_outgroup"]
ntrees = parsed_args["ntrees"]

# folders
script = relpath(abspath(@__FILE__), projectdir())
timestamp = now_string(minute=true)

parameters = @dict(
    potts_file, sample_potts_file,
    nleaves, nsweeps, b, outgroup, ntrees,
    script, timestamp,
)
identifier = savename(
    parsed_args["prefix"], parameters;
    ignores = ["potts_file", "sample_potts_file", "script", "timestamp"]
)
outfolder = joinpath(parsed_args["outfolder"], identifier)


@info "Using Potts model in $potts_file and eq. sample $sample_potts_file"
@info "Results saved in $outfolder"


dat_folder = joinpath(outfolder, "data")
if isdir(dat_folder)
    @info "Remove all contents of $dat_folder before simulating again, are you sure? [y/n]"
    yes = readline()
    !in(yes[1], ['Y', 'y']) && error("Aborting ($yes)")
    try
        rm(dat_folder; recursive=true)
    catch err
        @warn "Got error $err when trying to remove $dat_folder"
    end
end

open(joinpath(outfolder, "simulation_parameters.json"), "w") do f
    JSON3.pretty(f, JSON3.write(parameters))
end

## simulation
potts = DCAGraph(potts_file)
get_tree = () -> genealogy(YuleCoalescent(n, b))
ASRU.generate_trees(
    dat_folder, get_tree;
    add_outgroup=outgroup,
    outgroup_distance = :auto,
    M = ntrees,
)
foreach(ASRU.get_tree_folders(dat_folder)) do fol
    ASRU.simulate_sequences(fol, potts)
end
