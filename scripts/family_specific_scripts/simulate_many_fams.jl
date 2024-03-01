using DrWatson
@quickactivate "AutoRegressiveASR"

include(scriptsdir("simulate_data_potts_yule.jl"))
include(scriptsdir("asr_iqtree.jl"))
include(scriptsdir("asr_ardca.jl"))
include(scriptsdir("asr_ardca_diversity.jl"))
include(scriptsdir("asr_plm_at_root.jl"))


families = [
    "PF00076" => Dict(
        "prefix" => "PF00076",
        "potts" => datadir("Models", "PF00076/parameters_adaBM_PF00076.dat"),
        "sample_potts_eq" => datadir("Models", "PF00076/sample_adaBM_T200_PF00076.fasta"),
        "arnet" => datadir("Models", "PF00076/arnets/trained_on_potts/arnet_lJ0.001_lH0.0001.jld2"),
        "nsweeps" => 15,
    ),
    "PF00072" => Dict(
        "prefix" => "PF00072",
        "potts" => datadir("Models", "PF00072/parameters_dcatools_PF00072.dat"),
        "sample_potts_eq" => datadir("Models", "PF00072/sample_dcatools_T500_PF00072.fasta"),
        "arnet" => datadir("Models", "PF00072/arnets/trained_on_potts/arnet_lJ0.001_lH0.0001.jld2"),
        "nsweeps" => 15,
    ),
    "PF00014" => Dict(
        "prefix" => "PF00014",
        "potts" => datadir("Models", "PF00014/parameters_adaBM_PF00014.dat"),
        "sample_potts_eq" => datadir("Models", "PF00014/sample_adaBM_T200_PF00014.fasta"),
        "arnet" => datadir("Models", "PF00014/arnets/trained_on_potts/arnet_lJ0.001_lH0.0001.jld2"),
        "nsweeps" => 15,
    ),
    "PF00595" => Dict(
        "prefix" => "PF00595",
        "potts" => datadir("Models", "PF00595/parameters_adaBM_PF00595.dat"),
        "sample_potts_eq" => datadir("Models", "PF00595/sample_adaBM_T200_PF00595.fasta"),
        "arnet" => datadir("Models", "PF00595/arnets/trained_on_potts/arnet_lJ0.001_lH0.0001.jld2"),
        "nsweeps" => 15,
    ),
]

# Shared between families
nleaves = 50
ntrees = 25
add_outgroup = false
iqtree_model = "" # model finder if empty
outfolder = datadir("simulated", "potts_yule")


# to_simulate = ["PF00014"] # for testing
to_simulate = map(x -> x[1], families) # for all
for (name, fam) in filter(x -> x[1] in to_simulate, families)
    @info "Simulating for family $name"
    fam["nleaves"] = nleaves
    fam["ntrees"] = ntrees
    fam["add_outgroup"] = add_outgroup
    fam["outfolder"] = outfolder
    display(fam)

    folder = simulate_data_potts_yule(fam)
    fam["folder"] = folder

    fam["iqtree_model"] = iqtree_model
    asr_iqtree(fam)

    asr_ardca(fam)
    asr_ardca_sample_internals(fam)

    sample_internals_and_plm(fam)
end
