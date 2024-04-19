## simulate using ASR's simulator and models
function evolve(tree::Tree, model::ASR.EvolutionModel; kwargs...)
    ASR.Simulate.evolve(tree, model; kwargs...)
end

"""
    evolve(
        tree::Tree, graph::DCAGraph;
        alphabet = ASR.Alphabet(graph.mapping), leaves_fasta, internals_fasta, root_seq,
    )
"""
function evolve(
    tree::Tree, graph::DCAGraph;
    alphabet = ASR.Alphabet(graph.mapping),
    leaves_fasta="",
    internals_fasta="",
    root_seq = nothing,
    kwargs...
)
    @info "Simulating using `AutoRegressiveASR` functions"
    s0 = if isnothing(root_seq)
        DCATools.sample(graph, 1; Twait=1000) |> first |> vec
    else
        root_seq
    end
    leaf_sequences, internal_sequences = evolve_tree(
        tree, s0, graph, 1; kwargs...,
    )
    # write sequences to fasta if asked
    if !isempty(leaves_fasta)
        FASTAWriter(open(leaves_fasta, "w")) do writer
            for (name, seq) in leaf_sequences
                str_seq = ASR.intvec_to_sequence(seq; alphabet)
                write(writer, FASTARecord(name, str_seq))
            end
        end
    end
    if !isempty(internals_fasta)
        FASTAWriter(open(internals_fasta, "w")) do writer
            for (name, seq) in internal_sequences
                str_seq = ASR.intvec_to_sequence(seq; alphabet)
                write(writer, FASTARecord(name, str_seq))
            end
        end
    end

    return (leaf_sequences=leaf_sequences, internal_sequences = internal_sequences)
end

"""
    evolve(
        tree::Tree, graph::DCAGraph;
        alphabet = ASR.aa_alphabet, leaves_fasta, internals_fasta, root_seq,
    )
"""
function evolve(
    tree::Tree, arnet::ArDCA.ArNet;
    alphabet = ASR.aa_alphabet, leaves_fasta="", internals_fasta="", root_seq = nothing,
)
    @info "Simulating using `AutoRegressiveASR` functions -- ArNet model"
    s0 = isnothing(root_seq) ? ArDCA.sample(arnet, 1) |> vec : root_seq
    leaf_sequences, internal_sequences = evolve_tree(tree, s0, arnet, 1)
    # write sequences to fasta if asked
    if !isempty(leaves_fasta)
        FASTAWriter(open(leaves_fasta, "w")) do writer
            for (name, seq) in leaf_sequences
                str_seq = ASR.intvec_to_sequence(seq, alphabet)
                write(writer, FASTARecord(name, str_seq))
            end
        end
    end
    if !isempty(internals_fasta)
        FASTAWriter(open(internals_fasta, "w")) do writer
            for (name, seq) in internal_sequences
                str_seq = ASR.intvec_to_sequence(seq, alphabet)
                write(writer, FASTARecord(name, str_seq))
            end
        end
    end
end

#################################################################################
################################ Base functions #################################
#################################################################################

"""
    evolve_tree(tree::Tree, rootseq::AbstractVector{Int}, g, μ)

Return two dicts `leaf_seqs` and  `internal_seqs`. Does not modify `tree`.
"""
function evolve_tree(tree::Tree, rootseq::AbstractVector{Int}, g, μ; kwargs...)
    seqs = Dict(name => Int[] for name in map(label, nodes(tree)))
    _evolve!(tree.root, rootseq, g, μ, seqs; kwargs...)

    leaf_seqs = Dict(name => seqs[name] for name in map(label, leaves(tree)))
    internal_seqs = Dict(name => seqs[name] for name in map(label, internals(tree)))
    return leaf_seqs, internal_seqs
end

function _evolve!(node::TreeNode, seq::AbstractVector{Int}, g, μ, seqmap::Dict; kwargs...)
    seqmap[label(node)] = seq
    for c in children(node)
        c_seq = propagate(seq, μ*branch_length(c), g; kwargs...)
        _evolve!(c, c_seq, g, μ, seqmap)
    end
    return nothing
end
function _evolve!(
    node::TreeNode, seq::AbstractVector{Int}, g::DCAGraph, μ, seqmap::Dict; kwargs...
)
    seqmap[label(node)] = seq
    for c in children(node)
        τ = if !isinteger(μ*branch_length(c))
            τ = random_round(μ*branch_length(c))
            if abs(τ - μ*branch_length(c)) / τ > .1
                @warn "Rounding $(μ*branch_length(c)) to $(τ)"
            end
            τ
        else
            Int(μ*branch_length(c))
        end
        c_seq = propagate(seq, τ, g; kwargs...)
        _evolve!(c, c_seq, g, μ, seqmap)
    end
    return nothing
end

function random_round(x::Float64)
    f = Int(floor(x))
    dx = x - f
    if rand() > dx
        f
    else
        f + 1
    end
end
