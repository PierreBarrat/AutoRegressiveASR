using DrWatson
@quickactivate "AutoRegressiveASR"

using AutoRegressiveASR
using BackwardCoalescent
using Dates
using DCATools
using JSON3

function simulate_from_msa(parsed_args::AbstractDict; force=false)
    # msa of natural sequences
    in_fasta_file = parsed_args["msa_nat"]
    M = parsed_args[""]

    # teacher potts model
    potts_file = parsed_args["potts"] |> abspath
    sample_potts_file = parsed_args["sample_potts_eq"]
    potts = DCAGraph(potts_file)
end
