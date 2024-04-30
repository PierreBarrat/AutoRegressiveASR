@kwdef mutable struct ExpDist <: TreeTools.TreeNodeData
    depth::Float64 = 0.
    qdown::Union{Nothing, Float64} = nothing # info from below
    qup::Union{Nothing, Float64} = nothing # info from above
    Q::Float64 = 0.
end

"""
    weighted_leaf_distance!(r::TreeNode{ExpDist}, τ)
    weighted_leaf_distance!(t::Tree, τ)
"""
function weighted_leaf_distance!(r::TreeNode{ExpDist}, τ)
    up!(r, τ)
    down!(r, τ)
    depth!(r)
    return nothing
end
function weighted_leaf_distance!(t::Tree, τ)
    tree = convert(Tree{ExpDist}, t)
    weighted_leaf_distance!(tree.root, τ)
    return tree
end

"""
    set_to_zero!(n::TreeNode{ExpDist})
"""
function set_to_zero!(n::TreeNode{ExpDist})
    n.data.Q = 0.
    for c in n.child
        set_to_zero!(c)
    end
end

depth!(tree::Tree) = depth!(tree.root)
function depth!(n::TreeNode)
    if isleaf(n)
        n.data.depth = 0.
        return 0.
    end

    foreach(depth!, children(n))
    n.data.depth = let c = first(children(n))
        branch_length(c) + c.data.depth
    end
    return n.data.depth
end

function up!(n::TreeNode{ExpDist}, τ)
    if isleaf(n)
        n.data.qdown = 1.
        return 1.
    end

    n.data.qdown = sum(children(n)) do c
        t = branch_length(c)
        exp(-t/τ) * up!(c, τ)
    end
    n.data.qdown = min(1., n.data.qdown)
    return n.data.qdown
end

function down!(n::TreeNode{ExpDist}, τ)
    if isroot(n)
        n.data.qup = 0.
    else
        a = ancestor(n)
        n.data.qup = a.data.qup
        n.data.qup += sum(Iterators.filter(!=(n), children(a))) do m
            t = branch_length(m)
            exp(-t/τ) * m.data.qdown
        end
        n.data.qup *= exp(-branch_length(n) / τ)
        n.data.qup = min(1., n.data.qup)
    end
    foreach(c -> down!(c, τ), children(n))

    n.data.Q = min(1., n.data.qup + n.data.qdown)
    return n.data.qup
end

