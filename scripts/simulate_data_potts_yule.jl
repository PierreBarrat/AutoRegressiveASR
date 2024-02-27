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

using AutoRegressiveASR
using BackwardCoalescent
using DCATools
using JSON3

## Parameters
# Potts model
potts_file = parsed_args["potts"] |> abspath
sample_potts_file = parsed_args["sample_potts_eq"]
potts = DCAGraph(potts_file)

# Genealogy
nleaves = 50 # number of leaves
nsweeps = parsed_args["nsweeps"] # desired number of MCMC sweeps between root and leaves
L = potts.L
b = log(nleaves)/nsweeps/potts.L
outgroup = parsed_args["add_outgroup"]
ntrees = parsed_args["ntrees"]

# Out folder
timestamp = now_string(minute=true)

parameters = @dict(
    potts_file, sample_potts_file,
    nleaves, ntrees, nsweeps, b, outgroup,
    timestamp,
)
@tag!(parameters)
identifier = savename(
    parsed_args["prefix"], parameters;
    accesses = [:nleaves, :ntrees, :nsweeps, :outgroup],
    sort = true,
)
outfolder = joinpath(parsed_args["outfolder"], identifier)
mkpath(outfolder)

# Setting up folders
@info "Using Potts model in $(project_path(potts_file))"
@info "Results saved in $(project_path(outfolder))"

if isfile(sample_potts_file)
    cp(sample_potts_file, outfolder)
elseif isempty(sample_potts_file)
    @info "No sample file given"
else
    @warn "Could not find file $(abspath(sample_potts_file)) containing eq. sample of potts"
end

dat_folder = joinpath(outfolder, "data")
if isdir(dat_folder) || isfile(dat_folder)
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
get_tree = () -> genealogy(YuleCoalescent(nleaves, b))
ASRU.generate_trees(
    dat_folder, get_tree;
    add_outgroup=outgroup,
    outgroup_distance = :auto,
    M = ntrees,
)
foreach(ASRU.get_tree_folders(dat_folder)) do fol
    AutoRegressiveASR.simulate_sequences(fol, potts)
end
