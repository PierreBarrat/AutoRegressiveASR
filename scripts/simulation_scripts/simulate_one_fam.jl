using DrWatson
@quickactivate "AutoRegressiveASR"

include(scriptsdir("simulate_data_potts_yule.jl"))
include(scriptsdir("asr_iqtree.jl"))
include(scriptsdir("asr_ardca.jl"))
include(scriptsdir("asr_ardca_diversity.jl"))


fam_parameters = Dict(
    "prefix" => "PF00595",
    "potts" => datadir("Models", "PF00595/parameters_adaBM_PF00595.dat"),
    "sample_potts_eq" => datadir("Models", "PF00595/sample_adaBM_T200_PF00595.fasta"),
    "arnet" => datadir("Models", "PF00595/arnets/trained_on_potts/arnet_lJ0.001_lH0.0001.jld2"),

)
# Shared between runs
shared_parameters = Dict(
    "ntrees" => 25,
    "normalize_tree_height" => true,
    "add_outgroup" => false,
    "iqtree_model" => "PMB+I+G4", # model finder if empty
    "outfolder" => datadir("simulated", "potts_yule"),
)



# varying
parameters = Dict{Any, Any}(
    "nleaves" => [50],
    "nsweeps" => [2, 4, 8],
)


for p in dict_list(parameters)
    @info "Simulating for $p"
    # grouping parameters
    prm = convert(Dict{Any, Any}, p)
    merge!(prm, shared_parameters, fam_parameters)
    display(prm)

    folder = simulate_data_potts_yule(prm)
    prm["folder"] = folder

    asr_iqtree(prm)

    asr_ardca(prm)
    asr_ardca_sample_internals(prm)
end
