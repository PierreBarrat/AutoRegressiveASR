using DrWatson
quickactivate(@__DIR__, "AutoRegressiveASR")

using DelimitedFiles
using PlmDCA

get_fol(fam) = datadir("Models", fam)
families = Dict(
    "PF00014" => Dict("folder" => get_fol("PF00014"), "aln" => "PF00014_mgap6.fasta"),
    "PF00072" => Dict("folder" => get_fol("PF00072"), "aln" => "PF00072_mgap6_subsample.fasta"),
    "PF00076" => Dict("folder" => get_fol("PF00076"), "aln" => "PF00076_mgap6.fasta"),
    "PF00595" => Dict("folder" => get_fol("PF00595"), "aln" => "PF00595_mgap6.fasta"),
)

function write_scores(file, scores)
    scores = mapreduce(x -> collect(x)', vcat, scores) # to have it as i j score matrix
    open(file, "w") do io
        writedlm(io, vcat(["i" "j" "score"], scores), ' ')
    end
end

for (fam, dat) in families
    @info fam
    plmout = plmdca(joinpath(dat["folder"], dat["aln"]))
    write_scores(joinpath(dat["folder"], "scores_plm_on_nat.csv"), plmout.score)
end


