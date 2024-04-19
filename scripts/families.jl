using DrWatson
@quickactivate "AutoRegressiveASR"


families = Dict(
    "PF00076" => Dict(
        "prefix" => "PF00076",
        "L" => 70,
        "folder" => datadir("Models/PF00076/"),
        "potts" => datadir("Models", "PF00076/parameters_adaBM_PF00076.dat"),
        "sample_potts_eq" => datadir("Models", "PF00076/sample_adaBM_T200_PF00076.fasta"),
        "arnet" => datadir("Models", "PF00076/arnets/trained_on_potts/arnet_lJ0.01_lH0.001.jld2"),
        "sample_arnet_eq" => datadir("Models", "PF00076", "sample_arnet.fasta"),
        "aln_nat" => datadir("Models", "PF00076", "PF00076_mgap6.fasta"),
        "arnet_on_nat" => datadir("Models", "PF00076", "arnets/trained_on_nat/arnet_lJ0.01_lH0.001.jld2")

    ),
    "PF00072" => Dict(
        "prefix" => "PF00072",
        "L" => 112,
        "folder" => datadir("Models/PF00072/"),
        "potts" => datadir("Models", "PF00072/parameters_dcatools_PF00072.dat"),
        "sample_potts_eq" => datadir("Models", "PF00072/sample_dcatools_T500_PF00072.fasta"),
        "arnet" => datadir("Models", "PF00072/arnets/trained_on_potts/arnet_lJ0.01_lH0.001.jld2"),
        "sample_arnet_eq" => datadir("Models", "PF00072", "sample_arnet.fasta"),
        "aln_nat" => datadir("Models", "PF00072", "PF00072_mgap6_subsample.fasta"),
        "arnet_on_nat" => datadir("Models", "PF00072", "arnets/trained_on_nat/arnet_lJ0.01_lH0.001.jld2")
    ),
    "PF00014" => Dict(
        "prefix" => "PF00014",
        "L" => 53,
        "folder" => datadir("Models/PF00014/"),
        "potts" => datadir("Models", "PF00014/parameters_adaBM_PF00014.dat"),
        "sample_potts_eq" => datadir("Models", "PF00014/sample_adaBM_T200_PF00014.fasta"),
        "arnet" => datadir("Models", "PF00014/arnets/trained_on_potts/arnet_lJ0.01_lH0.001.jld2"),
        "sample_arnet_eq" => datadir("Models", "PF00014", "sample_arnet.fasta"),
        "aln_nat" => datadir("Models", "PF00014", "PF00014_mgap6.fasta"),
        "arnet_on_nat" => datadir("Models", "PF00014", "arnets/trained_on_nat/arnet_lJ0.01_lH0.001.jld2")
    ),
    "PF00595" => Dict(
        "prefix" => "PF00595",
        "L" => 82,
        "folder" => datadir("Models/PF00595/"),
        "potts" => datadir("Models", "PF00595/parameters_adaBM_PF00595.dat"),
        "sample_potts_eq" => datadir("Models", "PF00595/sample_adaBM_T200_PF00595.fasta"),
        "arnet" => datadir("Models", "PF00595/arnets/trained_on_potts/arnet_lJ0.01_lH0.001.jld2"),
        "sample_arnet_eq" => datadir("Models", "PF00595", "sample_arnet.fasta"),
        "aln_nat" => datadir("Models", "PF00595", "PF00595_mgap6.fasta"),
        "arnet_on_nat" => datadir("Models", "PF00595", "arnets/trained_on_nat/arnet_lJ0.01_lH0.001.jld2")
    ),
)
