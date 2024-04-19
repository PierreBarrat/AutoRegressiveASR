using DrWatson
@quickactivate("AutoRegressiveASR")

using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "folder"
        help = "Folder containing data -- itself containing `data/i/tree.nwk` etc..."
        required = true
    "arnet"
        help = "JLD2 file with ardca model (in the `:arnet` field)"
        required = true
end
parsed_args = parse_args(ARGS, s)

include(scriptsdir("asr_ardca.jl"))
asr_ardca(parsed_args)
