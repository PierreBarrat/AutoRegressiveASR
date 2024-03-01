using DrWatson
@quickactivate "AutoRegressiveASR"

using AutoRegressiveASR
using BackwardCoalescent
using Dates
using DCATools
using JSON3

function simulate_data_potts_yule(parsed_args::AbstractDict; force=false)
    ## Parameters
    # Potts model
    potts_file = parsed_args["potts"] |> abspath
    sample_potts_file = parsed_args["sample_potts_eq"]
    potts = DCAGraph(potts_file)

    # Genealogy
    nleaves = parsed_args["nleaves"] # number of leaves
    nsweeps = parsed_args["nsweeps"] # desired number of MCMC sweeps between root and leaves
    L = potts.L
    b = log(nleaves)/nsweeps/potts.L
    outgroup = parsed_args["add_outgroup"]
    ntrees = parsed_args["ntrees"]

    # Out folder
    timestamp = now()
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

    if isdir(outfolder) || isfile(outfolder)
        if force
            @info "Remove all contents of $outfolder before simulating again, are you sure? [y/n]"
            yes = readline()
            !occursin("yes", yes) && (@warn "Aborting simulation"; return outfolder)
            try
                rm(outfolder; recursive=true)
            catch err
                @warn "Got error $err when trying to remove $outfolder"
            end
        else
            @warn "$outfolder already exists, not simulating again"
            return outfolder
        end
    end
    mkpath(outfolder)

    # Setting up folders
    @info "Using Potts model in $(project_path(potts_file))"
    @info "Results saved in $(project_path(outfolder))"

    if isfile(sample_potts_file)
        try
            cp(sample_potts_file, joinpath(outfolder, basename(sample_potts_file)))
        catch err
            @info """Tried to copy $sample_potts_file to $outfolder but got error $err
            Maybe the file already existed"""
        end
    elseif isempty(sample_potts_file)
        @info "No sample file given"
    else
        @warn "Could not find file $(abspath(sample_potts_file)) containing eq. sample of potts"
    end


    open(joinpath(outfolder, "simulation_parameters.json"), "w") do f
        JSON3.pretty(f, JSON3.write(parameters))
    end

    ## simulation
    dat_folder = joinpath(outfolder, "data")
    get_tree = () -> genealogy(YuleCoalescent(nleaves, b))
    AutoRegressiveASR.generate_trees(
        dat_folder, get_tree;
        add_outgroup=outgroup,
        outgroup_distance = :auto,
        M = ntrees,
    )
    foreach(ASRU.get_tree_folders(dat_folder)) do fol
        AutoRegressiveASR.simulate_sequences(fol, potts)
    end

    @info "OUTPUT DIRECTORY $(abspath(outfolder))"
    return abspath(outfolder)
end
