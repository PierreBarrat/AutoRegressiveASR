"""
    reconstructed_to_df(
        tree_file, alignment_leaves, alignment_reconstructed;
        alignment_real = nothing,
        generative_model = nothing, # to compute loglikelihood of reconstructed sequences,
        model_consensus = nothing, # computed from `generative_model` if not provided
        inferred_tree_file = nothing,
        identifier, # identifier of the reconstruction, defaults to tree label,
        alphabet = :aa,
    )

Return a `DataFrame` with measurements about the quality of the reconstruction.
Input files are the tree used for ASR (or the real one) as a Newick file,
    the sequences at the leaves and the reconstructed ones at internal nodes,
    as `.fasta` files (expects `AbstractString`).

Optionally, the file containing real sequences at internal nodes (`alignment_real`) and the
    generative model used can be provided.
    For now, the latter should be and `ArDCA.ArNet`.
"""
function reconstructed_to_df(
    tree_file, alignment_leaves, alignment_reconstructed;
    alignment_real = nothing,
    generative_model = nothing,
    identifier = read_tree(tree_file).label,
    alphabet = :aa,
    inferred_tree_file = nothing, model_consensus = nothing, asr_state_file = nothing,
    include_reconstructed_sequence = false,
)
    #= Reading input data =#

    # tree
    tree = read_tree(tree_file)
    tree_inferred = isnothing(inferred_tree_file) ? nothing : read_tree(inferred_tree_file)

    # leaves
    leaf_sequences = FASTA.Reader(open(alignment_leaves, "r")) do reader
        map(rec -> description(rec) => sequence(rec), reader)
    end |> Dict

    # reconstructed at internals
    reconstructed = FASTA.Reader(open(alignment_reconstructed, "r")) do reader
        map(rec -> description(rec) => sequence(rec), reader)
    end |> Dict

    # real internals if provided
    real = if !isnothing(alignment_real)
        FASTA.Reader(open(alignment_real, "r")) do reader
            map(rec -> description(rec) => sequence(rec), reader)
        end |> Dict
    else
        nothing
    end

    # Checking we have the same internal nodes in all alignments and tree
    if !isnothing(real) && any(x -> !haskey(real, x), keys(reconstructed))
        error("""Alignment for real and reconstructed internals seems to differ.
            Got reconstructed: $(alignment_reconstructed) - and real: $(alignment_real)
        """)
    elseif any(x -> !haskey(reconstructed, x), internals(tree)) && any(x -> !in(x, tree), keys(reconstructed))
        @warn "Set of reconstructed sequences $(alignment_reconstructed) and tree internal nodes of $(tree_file) differ"
    end


    measures = MEASURES(generative_model)

    !include_reconstructed_sequence && delete!(measures, :reconstructed_sequence)

    # model consensus if needed
    model_cons = if haskey(measures, :hamming_to_model_consensus)
        if isnothing(model_consensus)
            consensus(generative_model; translate=true) # as a string
        else
            model_consensus
        end
    else
        nothing
    end
    aln_cons = consensus(leaf_sequences)

    # ground state if needed
    model_ground_state = if haskey(measures, :hamming_to_model_ground_state)
        ground_state(generative_model; translate=true) # as a string
    else
        nothing
    end

    # state file for entropy
    asr_state_profiles = if !isnothing(asr_state_file)
        parse_asr_state_file(asr_state_file)
    else
        nothing
    end

    # Arranging all data in one tuple for convenient call of functions
    data = (
        tree = tree,
        tree_inferred = tree_inferred,
        leaf_sequences = leaf_sequences,
        reconstructed_sequences = reconstructed,
        real_sequences = real,
        model = generative_model,
        alphabet = alphabet,
        model_consensus = model_cons,
        model_ground_state = model_ground_state,
        alignment_consensus = aln_cons,
        asr_state_profiles = asr_state_profiles,
    )

    #= Measures =#
    labels = @chain begin
        Iterators.filter(isinternal, POT(tree))
        map(label, _)
        filter!(x -> haskey(reconstructed, x), _) # only consider nodes that have been reconstructed
    end
    df = DataFrame(
        identifier = map(_ -> identifier, labels),
        label = labels,
        sequence = map(lab -> data.reconstructed_sequences[lab], labels),
    )
    for (name, measure) in measures
        df[!, name] = map(lab -> measure(data, lab), labels)
    end
    select!(df, Cols(Not(:sequence), :sequence)) # change order of cols for nicer display
    sort!(df, [:identifier, :label])
    return  df
end


#=================================================================================#
################################### Measures ######################################
#=================================================================================#

#=
All of these functions take `data` and a single `label` as input, and return a number.
=#


node_depth(data, label) = TreeTools.distance_to_closest_leaf(data.tree, label)
function node_depth_inferred(data, label)
    return if isnothing(data.tree_inferred)
        missing
    else
        TreeTools.distance_to_closest_leaf(data.tree_inferred, label)
    end
end

reconstructed_sequence(data, label) = data.reconstructed_sequences[label]

myloglikelihood(S, g::DCAGraph) = -DCATools.energy(g, S)
myloglikelihood(S, arnet::ArNet) = ArDCA.loglikelihood(S, arnet)
function loglikelihood(data, label)
    if isnothing(data.model)
        return missing
    end

    S = ASR.sequence_to_intvec(
        data.reconstructed_sequences[label], alphabet = data.alphabet
    )
    return myloglikelihood(S, data.model)
end



function hamming_to_closest_leaf(data, lab)
    node = data.tree[lab]
    leaf = argmin(Iterators.map(lf -> label(lf) => lf, leaves(data.tree))) do (_, leaf)
        distance(node, leaf)
    end[1]
    hamming(data.leaf_sequences[leaf], data.reconstructed_sequences[lab])
end

function hamming_to_root(data, lab)
    return if isnothing(data.real_sequences) ||
        !haskey(data.reconstructed_sequences, label(data.tree.root))
        missing
    else
        hamming(
            data.reconstructed_sequences[label(data.tree.root)],
            data.reconstructed_sequences[lab]
        )
    end
end

function hamming_to_real(data, label)
    return if isnothing(data.real_sequences)
        missing
    else
        hamming(data.real_sequences[label], data.reconstructed_sequences[label])
    end
end

function hamming_to_real_nogap(data, label)
    return if isnothing(data.real_sequences)
        missing
    else
        hamming(
            data.real_sequences[label], data.reconstructed_sequences[label];
            exclude_gaps=true,
        )
    end
end

function hamming_to_model_consensus(data, label)
    return if isnothing(data.model_consensus)
        missing
    else
        hamming(data.model_consensus, data.reconstructed_sequences[label])
    end
end
function hamming_to_model_ground_state(data, label)
    return if isnothing(data.model_ground_state)
        missing
    else
        hamming(data.model_ground_state, data.reconstructed_sequences[label])
    end
end

function hamming_to_aln_consensus(data, label)
    hamming(data.alignment_consensus, data.reconstructed_sequences[label])
end

function entropy(data, label)
    return if isnothing(data.asr_state_profiles)
        missing
    else
        ASR.entropy(data.asr_state_profiles[label].model)
    end

end

function _isroot(data, label)
    return data.tree[label].isroot
end

MEASURES(generative_model) = Dict(
    :node_depth => node_depth,
    :node_depth_inferred => node_depth_inferred,
    :root => _isroot,
    :loglikelihood => loglikelihood,
    :hamming_to_real => hamming_to_real,
    :hamming_to_real_nogap => hamming_to_real_nogap,
    :hamming_to_closest_leaf => hamming_to_closest_leaf,
    :hamming_to_root => hamming_to_root,
    :hamming_to_model_consensus => hamming_to_model_consensus,
    :hamming_to_model_ground_state => hamming_to_model_ground_state,
    :hamming_to_aln_consensus => hamming_to_aln_consensus,
    :entropy => entropy,
    :reconstructed_sequence => reconstructed_sequence,
)
MEASURES(::DCAGraph) = Dict(
    :node_depth => node_depth,
    :node_depth_inferred => node_depth_inferred,
    :root => _isroot,
    :loglikelihood => loglikelihood,
    :hamming_to_real => hamming_to_real,
    :hamming_to_real_nogap => hamming_to_real_nogap,
    :hamming_to_closest_leaf => hamming_to_closest_leaf,
    :hamming_to_root => hamming_to_root,
    :hamming_to_model_consensus => hamming_to_model_consensus,
    # :hamming_to_model_ground_state => hamming_to_model_ground_state,
    :hamming_to_aln_consensus => hamming_to_aln_consensus,
    :entropy => entropy,
    :reconstructed_sequence => reconstructed_sequence,
)
MEASURES(::Nothing) = Dict(
    :node_depth => node_depth,
    :node_depth_inferred => node_depth_inferred,
    :root => _isroot,
    # :loglikelihood => loglikelihood,
    :hamming_to_real => hamming_to_real,
    :hamming_to_real_nogap => hamming_to_real_nogap,
    :hamming_to_closest_leaf => hamming_to_closest_leaf,
    :hamming_to_root => hamming_to_root,
    # :hamming_to_model_consensus => hamming_to_model_consensus,
    # :hamming_to_model_ground_state => hamming_to_model_ground_state,
    :hamming_to_aln_consensus => hamming_to_aln_consensus,
    :reconstructed_sequence => reconstructed_sequence,
)
