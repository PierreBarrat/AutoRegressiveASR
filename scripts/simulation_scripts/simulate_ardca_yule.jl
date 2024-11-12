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
    "iqtree_model" => "best", # model finder if empty
    "opt_bl" => :scale,
    "outfolder" => datadir("simulated", "arnet_yule"),
    "remove_iqtree_statefile" => false, # to not use up too much space -- set to false if diversity needed
    "prefix" => "",
    "suffix" => "scaleiqtree_clean",
)

# varying
parameters = Dict{Any, Any}(
    "nleaves" => [100],
    "treeheight" => [2.],
)


# to_simulate = ["PF00014"]
to_simulate = collect(keys(families))
for (fam, fam_parameters) in families
    !in(fam, to_simulate) && continue
    @info "Simulating $fam"

    for p in dict_list(parameters)
        @info "Simulating for $p"
        # grouping parameters
        prm = convert(Dict{Any, Any}, p)
        merge!(prm, shared_parameters)
        if prm["iqtree_model"] == "best"
            prm["iqtree_model"] = fam_parameters["iqtree_model_ardca"]
        end
        prm["prefix"] = fam_parameters["prefix"]
        prm["generative_model"] = fam_parameters["arnet"]
        prm["sample_equilibrium"] = fam_parameters["sample_arnet_eq"]

        display(prm)

        folder = simulate_data_ardca_yule(prm)
        prm["folder"] = folder

        #= With base iqtree model =#
        _, iqtree_prefix = asr_iqtree(prm)
        #= With iqtree profile models =#
        ref_iqtree_model = prm["iqtree_model"]
        prm["iqtree_model"] = ref_iqtree_model * "+C10"
        asr_iqtree(prm)
        prm["iqtree_model"] = ref_iqtree_model * "+C60"
        asr_iqtree(prm)

        prm["arnet"] = fam_parameters["arnet"] # for reconstruction
        prm["iqtree_prefix"] = iqtree_prefix # for scaling branch lengths
        asr_ardca(prm)
        if prm["treeheight"] == 2.
            asr_ardca_sample_internals(prm)
        end
    end
end
