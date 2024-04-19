using DrWatson
@quickactivate "AutoRegressiveASR"

include(scriptsdir("simulate_data_potts_yule.jl"))
include(scriptsdir("asr_iqtree.jl"))
include(scriptsdir("asr_ardca.jl"))
include(scriptsdir("asr_ardca_diversity.jl"))

include(scriptsdir("families.jl"))

# Shared between runs
shared_parameters = Dict(
    "ntrees" => 100,
    "nsim_per_tree" => 1,
    "normalize_tree_height" => true,
    "add_outgroup" => false,
    "iqtree_model" => "", # model finder if empty
    "asr_opt_bl" => :opt,
    "outfolder" => datadir("simulated", "potts_yule"),
)

# varying
parameters = Dict{Any, Any}(
    "nleaves" => [100],
    "nsweeps" => [8],
)


to_simulate = ["PF00076", "PF00595"]
for (fam, fam_parameters) in families
    !in(fam, to_simulate) && continue
    @info "Simulating $fam"
    for p in dict_list(parameters)
        @info "Simulating for $p"
        # grouping parameters
        prm = convert(Dict{Any, Any}, p)
        merge!(prm, shared_parameters)
        prm["prefix"] = fam_parameters["prefix"]
        prm["generative_model"] = fam_parameters["potts"]
        prm["sample_equilibrium"] = fam_parameters["sample_potts_eq"]
        display(prm)

        folder = simulate_data_potts_yule(prm)
        prm["folder"] = folder
        asr_iqtree(prm)

        prm["arnet"] = fam_parameters["arnet"] # for reconstruction -- always arnet
        asr_ardca(prm)
        asr_ardca_sample_internals(prm)
    end
end
