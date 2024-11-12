using DrWatson
@quickactivate "AutoRegressiveASR"

using CSV
using DataFrames
using Glob

include(scriptsdir("families.jl"))
includet(scriptsdir("root_position/simulate_arnet.jl"))
includet(scriptsdir("root_position/reconstruct_ardca.jl"))
includet(scriptsdir("root_position/gather_results.jl"))

# Shared between runs
parameters = Dict{String,Any}(
    "ntrees" => 10,
    "nleaves" => 100,
    "treeheight" => 2.,
    "add_outgroup" => false,
    "opt_bl" => :real,
    "outfolder" => datadir("simulated", "arnet_yule_rootposition"),
    "prefix" => "",
    "suffix" => "",
    "family" => "PF00072",
)

parameters["prefix"] = parameters["family"]
parameters["generative_model"] = families[parameters["family"]]["arnet"]

parameters["famfolder"] = simulate_data_arnet(parameters) # outfolder/PFXXX...
# Simulating with ardca - measure things - store in dataframe
for tree_folder in glob(["simulations", r"[0-9]+"], parameters["famfolder"])
    @info tree_folder
    parameters["dat_folder"] = joinpath(tree_folder, "data")
    asr_ardca(parameters)
    data = gather_results(parameters["dat_folder"])
    CSV.write(joinpath(parameters["dat_folder"], "../measures.csv"), data)
end
# group the dataframes in a big one
df = mapreduce(vcat, glob(["simulations", r"[0-9]+"], parameters["famfolder"])) do fol
    CSV.read(joinpath(fol, "measures.csv"), DataFrame)
end
CSV.write(joinpath(parameters["famfolder"], "measures.csv"), df)






