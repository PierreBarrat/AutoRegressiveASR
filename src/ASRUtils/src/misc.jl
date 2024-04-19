"""
    get_tree_folders(dir)

Find and return all directories of the form "dir/i" where `i` is an integer.
"""
function get_tree_folders(dir)
    return filter(readdir(dir; join=true)) do f
        x = tryparse(Int64, basename(f))
        isnothing(x) ? false : true
    end
end

function hamming(X::AbstractString, Y::AbstractString; normalize=true, exclude_gaps = false)
    length(X) != length(Y) && error("""Length differ
        Got $X of length $(length(X)) and $Y of length $(length(Y))
    """)
    if !exclude_gaps
        S = sum(zip(X,Y)) do (x,y) x!= y end
        return normalize ? S / length(X) : S
    end

    S, Z = sum(zip(X, Y)) do (x, y)
        r = [0, 0]
        if x != '-' && y != '-'
            r[2] += 1
            x != y && (r[1] += 1)
        end
        r
    end
    return normalize ? S/Z : S
end
function hamming(
    X::AbstractVector, Y::AbstractVector;
    normalize = true, exclude_gaps = false, gap_state::Union{Int, Nothing} = nothing
)
    length(X) != length(Y) && error("""Length differ
        Got $X of length $(length(X)) and $Y of length $(length(Y))
    """)
    if !exclude_gaps
        S = sum(zip(X,Y)) do (x,y) x!= y end
        return normalize ? S / length(X) : S
    end

    isnothing(gap_state) && error("Provide an integer `gap_state` as kwarg")
    S, Z = sum(zip(X, Y)) do (x, y)
        r = [0, 0]
        if x != gap_state && y != gap_state
            r[2] += 1
            x != y && (r[1] += 1)
        end
        r
    end
    return normalize ? S/Z : S
end

"""
    balanced_tree(K, T)

Generate a balanced bifurcating with `K` "floors" (*i.e.* 2^K leaves) and time to root `T`.
"""
function balanced_tree(K, T; branch_length = :constant, randeps = 0.)
    n = 2^K

    t(k) = if branch_length == :constant
        T/K * (1+randn()*randeps)
    elseif branch_length == :exponential
        T*1/2^k * (1+randn()*randeps)
    elseif branch_length == :linear
        2*(K-k+1)*T / (K*(K+1)) * (1+randn()*randeps)
    else
        error("Unknown value for `branch_length` kwarg")
    end

    tree = basic_shape(t(1))
    for k in 2:K
        for l in collect(leaves(tree))
            sub = basic_shape(t(k))
            graft!(tree, sub.root, l; graft_on_leaf=true)
            branch_length!(tree[label(sub.root)], 0)
        end
        TreeTools.remove_internal_singletons!(tree)
    end
    return tree
end

function ladder_tree(K, T)
    dt = T/K
    tree = ASRU.basic_shape(dt)
    left_most = first(children(tree.root))

    for k in 2:K
        # increase all branch lengths but the left most one
        foreach(Iterators.filter(!=(left_most), leaves(tree))) do node
            branch_length!(node, branch_length(node) + dt)
        end

        sub = ASRU.basic_shape(dt)
        graft!(tree, sub.root, left_most; graft_on_leaf=true)
        branch_length!(tree[label(sub.root)], 0)
        left_most = first(children(sub.root))

        TreeTools.remove_internal_singletons!(tree)
    end

    tree
end

# helper function for balanced_tree
function basic_shape(t::Number)
    l1, l2, r = [TreeTools.make_random_label() for _ in 1:3]
    return parse_newick_string("($l1:$t,$l2:$t)$r;")
end


function consensus(arnet::ArNet; alphabet = ASR.aa_alphabet, translate=true, kwargs...)
    return @chain begin
        ASR.ProfileModel(arnet; kwargs...)
        ASR.ml_sequence
        translate ? ASR.intvec_to_sequence(_; alphabet) : _
    end
end

function consensus(
    graph::DCAGraph;
    alphabet = ASR.Alphabet(graph.mapping), translate=true, M = 1000, kwargs...
)
    return @chain graph begin
        DCATools.sample(_, M; Twait = 25)
        DCATools.consensus
        first
        translate ? ASR.intvec_to_sequence(_; alphabet) : _
    end
end

function consensus(
    tree::Tree{ASR.AState{q}}; alphabet = ASR.default_alphabet(q), translate=true
) where q
    leaf_sequences = map(x -> x.data.sequence, leaves(tree))
    L = length(first(leaf_sequences))
    cons = map(1:L) do i
        [s[i] for s in leaf_sequences] |> countmap |> argmax
    end
    return translate ? ASR.intvec_to_sequence(cons; alphabet) : cons
end

function consensus(leaf_sequences::AbstractDict)
    L = length(first(values(leaf_sequences)))
    cons = map(1:L) do i
        [s[i] for s in values(leaf_sequences)] |> countmap |> argmax
    end
    toseq(x::Vector{Int}) = x
    toseq(x) = prod(x) # if vector of char or string
    return toseq(cons)
end

"""
    easy_smooth(x, y; w = 10)
    easy_smooth(df::DataFrame, f1, f2; kwargs...)
"""
function easy_smooth(x, y; w = 25, alg=:sma, outliers_left = 0., outliers_right = .05)
    idx = exclude_outliers_perm(x; left=outliers_left, right=outliers_right)
    xs = x[idx]
    ys = y[idx]
    (xsmth, ysmth, ystd, nsamples) = if alg == :loess
        xs, loess(xs, ys; q=w)(xs), nothing, nothing
    elseif alg == :sma
        xs, sma(ys, w, true), nothing, nothing
    elseif alg == :hist
        hist_smooth(xs, ys; nbins=w, outliers_left=0., outliers_right=0.) #outliers already removed
    else
        error("Unknwon smoothing alg $alg. ")
    end
    return collect.(skipmissings(xsmth, ysmth, ystd, nsamples))
    # return collect.(skipmissings(xs, sma(ys, w, true)))
end
easy_smooth(df::DataFrame, f1, f2; kwargs...) = easy_smooth(df[!, f1], df[!, f2]; kwargs...)

function hist_smooth(x, y; nbins=100, outliers_left=0., outliers_right = .05)
    idx = exclude_outliers_perm(x; left=outliers_left, right=outliers_right)
    xs = x[idx]
    ys = y[idx]

    xmin, xmax = extrema(xs)
    edges = range(xmin, xmax, length=nbins+1)
    summed_weights = fit(Histogram, xs, edges).weights |> cumsum

    y_av = map(1:nbins) do i
        ilow = i == 1 ? 1 : summed_weights[i-1]
        mean(ys[ilow:summed_weights[i]])
    end
    y_std = map(1:nbins) do i
        ilow = i == 1 ? 1 : summed_weights[i-1]
        std(ys[ilow:summed_weights[i]])
    end
    nsamples = map(1:nbins) do i
        ilow = i == 1 ? 1 : summed_weights[i-1]
        length(ys[ilow:summed_weights[i]])
    end

    bin_centers = (edges[1:end-1] .+ edges[2:end])/2
    return bin_centers, y_av, y_std, nsamples
end

function exclude_outliers_perm(x; left = 0., right = 0.)
    idx = sortperm(x)
    L = length(x)
    rb = L - Int(floor(right * L))
    lb = 1+Int(floor(left * L))
    return idx[lb:rb]
end

function ground_state(arnet::ArNet; alphabet = ASR.aa_alphabet, translate = true)
    # @extract arnet:H J p0 idxperm
    q = length(arnet.p0)
    N = length(arnet.H) # here N is N-1 !!
    backorder = sortperm(arnet.idxperm)

    totH = Vector{Float64}(undef, q)
    ground_state = Vector{Int}(undef, N + 1)
    ground_state[1] = argmax(arnet.p0)
    for site in 1:N
        Js = arnet.J[site]
        h = arnet.H[site]
        copy!(totH, h)
        for i in 1:site
            for a in 1:q
                totH[a] += Js[a, ground_state[i], i]
            end
        end
        ground_state[site+1] = argmax(totH)
    end

    ArDCA.permuterow!(ground_state, backorder)

    return if translate
        ASR.intvec_to_sequence(ground_state; alphabet)
    else
        ground_state
    end
end

"""
    rearrange_from_model(tree, model)

Assume `tree` and `model` have the same topology when unrooted. This will
- reroot `tree` in the same way as `model`,
- relabel internal nodes of `tree` in the same way as model
Both trees need to *share leaf node labels* in order for this to work.
Return a copy of `tree`.
"""
function rearrange_from_model(input_tree, model)
    @assert share_labels(input_tree, model) "Only works on trees that share leaf nodes"

    tree = copy(input_tree)
    TreeTools.root!(tree; method=:model, model)

    @assert SplitList(tree) == SplitList(model) "Trees differ in topology, aborting."
    S = SplitList(tree)
    for i in 1:length(S)
        nmodel = model[S, i]
        n = tree[S, i]
        label(nmodel) != label(n) && label!(tree, n, label(nmodel))
    end
    return tree
end
