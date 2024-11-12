using DrWatson
@quickactivate "AutoRegressiveASR"

using AutoRegressiveASR
using Chain
using CSV
using DataFrames
using FASTX
using Glob
using TreeTools


function gather_results(folder)
    real_internals = fasta_to_dict(joinpath(folder, "alignment_internals.fasta"))
    original_reconstructed_internals = fasta_to_dict(
        joinpath(folder, "original_reconstructed_internals_ML.fasta")
    )
    original_root = @chain folder joinpath("tree.nwk") read_tree root label

    mapreduce(vcat, glob([r"[0-9]+"], folder)) do subfol
        results(subfol, real_internals, original_reconstructed_internals, original_root)
    end
end

function results(
    folder,
    real_internals::AbstractDict,
    original_reconstructed_internals::AbstractDict,
    original_root::AbstractString
)
    tree = read_tree(joinpath(folder, "tree.nwk"); check=false) # warning because root is singleton --> check=false
    reconstructed_internals = fasta_to_dict(
        joinpath(folder, "reconstructed_internals_ML.fasta")
    )
    # values below are the same for all nodes
    folder_label = hash(folder)
    current_root = label(root(tree))
    delta_root = distance(tree, current_root, original_root) # how re-rooted this tree is
    parameters = Dict(
        "reconstructed_internals" => reconstructed_internals,
        "real_internals" => real_internals,
        "original_reconstructed_internals" => original_reconstructed_internals,
    )
    # node measures
    node_iter = traversal(tree, :postorder; leaves=false) do node
        label(node) != original_root # original root is a singleton in other trees --> it will have been removed
    end

    df = DataFrame()
    for (name, measure) in MEASURES
        node_iter = traversal(tree, :postorder; leaves=false) do node
            label(node) != original_root # original root is a singleton in other trees --> it will have been removed
        end
        df[!, name] = map(node -> measure(label(node), parameters), node_iter)
    end

    n_internals = size(df, 1)
    df[!, :identifier] = [folder_label for _ in 1:n_internals]
    df[!, :current_root] = [current_root for _ in 1:n_internals]
    df[!, :delta_root] = repeat([delta_root], n_internals)

    # sort for nicer display
    select!(df, Cols(:identifier, :node, :delta_root, Not(:sequence), :sequence)) # change order of cols for nicer display
    sort!(df, [:identifier, :node])

    return df
end


get_node_label(node, _) = node # already passing label
get_sequence(node, parameters) = parameters["reconstructed_internals"][node]

function hamming_to_real(node, parameters)
    @unpack reconstructed_internals, real_internals = parameters
    return ASRU.hamming(reconstructed_internals[node], real_internals[node])
end
function hamming_to_original_reconstruction(node, parameters)
    @unpack reconstructed_internals, original_reconstructed_internals = parameters
    return ASRU.hamming(reconstructed_internals[node], original_reconstructed_internals[node])
end
function original_hamming_to_real(node, parameters)
    @unpack original_reconstructed_internals, real_internals = parameters
    return ASRU.hamming(original_reconstructed_internals[node], real_internals[node])
end

function fasta_to_dict(fastafile)
    return FASTA.Reader(open(fastafile, "r")) do reader
        map(rec -> description(rec) => sequence(rec), reader)
    end |> Dict
end

# Only measures that need the node in here
MEASURES = Dict(
    :node => get_node_label,
    :sequence => get_sequence,
    :hamming_to_real => hamming_to_real,
    :hamming_to_original_reconstruction => hamming_to_original_reconstruction,
    :original_hamming_to_real => original_hamming_to_real,
)
