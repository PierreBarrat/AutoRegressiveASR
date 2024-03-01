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
    "--prefix"
        help = "prefix for output folder"
        arg_type = String
        default = ""
    "--add_outgroup"
        help = "Add outgroup for rooting"
        action = :store_true
    "--nsweeps"
        help = "desired number of MCMC sweeps between root and leaves"
        arg_type = Int
        default = 1
    "--ntrees"
        help = "Number of trees to simulate"
        arg_type = Int
        default = 1
    "--nleaves"
        help = "Number of leaves in the trees"
        arg_type = Int
        default = 50
    "potts"
        help = "file containing the potts model"
        arg_type = String
        required = true
end
parsed_args = parse_args(ARGS, s)

include(scriptsdir("simulate_data_potts_yule.jl"))

simulate_data_potts_yule(parsed_args)
