using DrWatson
@quickactivate "AutoRegressiveASR"

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--model"
        help = "Model for iqtree, *e.g.* `Blosum62+I+G4`. If empty, use model finder."
        default = ""
    "folder"
        help = "Folder containing data, with subfolder `data/i/alignment_leaves.fasta` and `data/i/tree.nwk`"
        arg_type = String
        required = true
end
parsed_args = parse_args(ARGS, s)

include(scriptsdir("asr_iqtree.jl"))

asr_iqtree(parsed_args)
