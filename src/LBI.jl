"""
    mutable struct LBIData <: TreeTools.TreeNodeData

Data used to compute the Local Branching Index.
"""
mutable struct LBIData <: TreeTools.TreeNodeData
    message_down::Float64
    message_up::Float64
    lbi::Float64
end
function LBIData(;
    message_down=0.,
    message_up=0.,
    LBI=0.,
)
    return LBIData(message_down, message_up, LBI)
end


let state::Float64 = 0.
    global lbi_newstate(x) = (state=x)
    global lbi_getstate() = state
    global lbi_resetstate() = lbi_newstate(0)
end

"""
    lbi!(r::TreeNode{LBIData}, τ)
    lbi!(t::Tree{LBIData}, τ)
"""
function lbi!(r::TreeNode{LBIData}, τ::Real ; normalize=false)
    get_message_up!(r::TreeNode{LBIData}, τ)
    send_message_down!(r::TreeNode{LBIData}, τ)
    compute_lbi!(r::TreeNode{LBIData})
    if normalize
        normalize_lbi!(r)
    end
    return nothing
end
function lbi!(t::Tree, τ; normalize=false)
    tree = convert(Tree{LBIData}, t)
    lbi!(tree.root, τ, normalize=normalize)
    return tree
end


"""
    set_to_zero!(n::TreeNode{LBIData})
"""
function set_to_zero!(n::TreeNode{LBIData})
    n.data.lbi = 0.
    for c in n.child
        set_to_zero!(c)
    end
end

"""
    compute_lbi!(n::TreeNode{LBIData})
"""
function compute_lbi!(n::TreeNode{LBIData})
    n.data.lbi = n.data.message_down
    for c in n.child
        compute_lbi!(c)
        n.data.lbi += c.data.message_up
    end
end

"""
    get_message_up!(n::TreeNode{LBIData}, τ)

Get message going up from node `n`. Ask for up messages `m_c` for all children `c` of `n`. Return `exp(-t/τ) * sum(m_c) + τ * (1 - exp(-t/τ)` where `t=n.tau`. Field `message_up` in `n`'s data is modified.
"""
function get_message_up!(n::TreeNode{LBIData}, τ)
    n.data.message_up = 0.
    n.data.message_down = 0.
    for c in n.child
        n.data.message_up += get_message_up!(c, τ)
    end
    if !isroot(n)
        t = branch_length(n)
        n.data.message_up *= exp(-t/τ)
        n.data.message_up += τ*(1 - exp(-t/τ))
    end
    return n.data.message_up
end

"""
    send_message_down!(n::TreeNode{LBIData}, τ)

Send message going down from node `n`. Field `c.message_down` is modified for all children `c` of `n`. No return value.

## Note
It's a bit easier to think of this function as operating on `c.anc` where `c` is a child of `n`. Maybe I should code it like this.
"""
function send_message_down!(n::TreeNode{LBIData}, τ)
    for c1 in children(n)
        c1.data.message_down = n.data.message_down
        for c2 in Iterators.filter(!=(c1), children(n))
            c1.data.message_down += c2.data.message_up
        end
        t = branch_length(c1)
        c1.data.message_down *= exp(-t/τ)
        c1.data.message_down += τ*(1 - exp(-t/τ))

        send_message_down!(c1,τ)
    end
    return nothing
end

"""
    normalize_lbi!(r::TreeNode{LBIData})
    normalize_lbi!(t::Tree)
"""
function normalize_lbi!(r::TreeNode{LBIData})
    max_lbi = get_max_lbi(r::TreeNode{LBIData})
    normalize_lbi!(r::TreeNode{LBIData}, max_lbi)
end
normalize_lbi!(t::Tree) = normalize_lbi!(t.root)
function normalize_lbi!(n::TreeNode{LBIData}, max_lbi::Float64)
    if max_lbi == 0.
        @error "Maximum LBI is 0, cannot normalize."
    end
    for c in n.child
        normalize_lbi!(c, max_lbi)
    end
    n.data.lbi /= max_lbi
end

"""
    get_max_lbi(r::TreeNode{LBIData})

Return maximal LBI value of the subtree below `r`.
"""
function get_max_lbi(r::TreeNode{LBIData})
    lbi_resetstate()
    _get_max_lbi(r)
    return lbi_getstate()
end
function _get_max_lbi(n::TreeNode{LBIData})
    for c in n.child
        _get_max_lbi(c)
    end
    if n.data.lbi > lbi_getstate()
        lbi_newstate(n.data.lbi)
    end
end

