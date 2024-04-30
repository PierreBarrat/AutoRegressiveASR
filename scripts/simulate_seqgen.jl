using BioSequenceMappings
using Dates
using DelimitedFiles
using JLD2
using Random
using StatsBase
using TreeTools

include(scriptsdir("families.jl"))


function line_tree(tvals)
    tree = Tree(TreeNode(label="init"))
    for t in filter(!=(0), tvals)
        n = TreeNode(label="t$(t)")
        graft!(tree, n, first(leaves(tree)); tau = t, graft_on_leaf = true)
    end
    return tree
end
