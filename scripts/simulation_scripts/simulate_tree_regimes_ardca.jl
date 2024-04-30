using DrWatson
@quickactivate "AutoRegressiveASR"

include(scriptsdir("generate_data_ardca_coalescent.jl"))
include(scriptsdir("asr_iqtree.jl"))
include(scriptsdir("asr_ardca.jl"))

include(scriptsdir("families.jl"))


# Shared between runs
shared_parameters = Dict(
    "ntrees" => 100,
    "nsim_per_tree" => 1,
    "normalize_tree_height" => true,
    "add_outgroup" => false,
    "iqtree_model" => "PMB+I+G4", # model finder if empty
    "asr_opt_bl" => :fromreal,
    "outfolder" => datadir("simulated", "tree_regimes_arnet"),
    "remove_iqtree_statefile" => true # to not use up too much space -- set to false if diversity needed
)

# varying
parameters = Dict{Any, Any}(
    "nleaves" => [50],
    "treeheight" => [.25, .5, .75, 1., 1.5, 2., 3., 4.,],
    "coalescent" => ["Kingman", "Yule"]
)


to_simulate = ["PF00595"]
for (fam, fam_parameters) in families
    !in(fam, to_simulate) && continue
    @info "Simulating $fam"
    for p in dict_list(parameters)
        @info "Simulating for $p"
        # grouping parameters
        prm = convert(Dict{Any, Any}, p)
        merge!(prm, shared_parameters)
        prm["prefix"] = fam_parameters["prefix"]
        prm["generative_model"] = fam_parameters["arnet"]
        prm["sample_equilibrium"] = fam_parameters["sample_arnet_eq"]

        display(prm)

        folder = simulate_data_ardca_coalescent(prm)
        prm["folder"] = folder

        asr_iqtree(prm)

        prm["arnet"] = fam_parameters["arnet"] # for reconstruction -- always arnet
        asr_ardca(prm)
    end
end
