using DrWatson
@quickactivate "AutoRegressiveASR"

include(scriptsdir("generate_data_ardca_yule.jl"))
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
    "asr_opt_bl" => :fromreal,
    "outfolder" => datadir("simulated", "arnet_yule"),
    "remove_iqtree_statefile" => false # to not use up too much space -- set to false if diversity needed
)

# varying
parameters = Dict{Any, Any}(
    "nleaves" => [100],
    "treeheight" => [2.],
)


# to_simulate = ["PF00014"]#, "PF00595"]
to_simulate = filter(!=("PF00014"), collect(keys(families)))
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

        folder = simulate_data_ardca_yule(prm)
        prm["folder"] = folder

        asr_iqtree(prm)

        prm["arnet"] = fam_parameters["arnet"] # for reconstruction -- always arnet
        asr_ardca(prm)
        # if prm["treeheight"] == 2.
        asr_ardca_sample_internals(prm)
        # end
    end
end
