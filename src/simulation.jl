"""
    generate_trees(outdir; n=50, b=1, M=10)

For `r` in `1:M` generate a hierarchy of folders:
```
outdir/r/tree.nwk
```
with `tree` sampled from the Yule coalescent with parameter `b` and `n` lineages.
The average height of trees obtained like this is `log(n)/b`.
The number of internal nodes distant of `x` from the leaves goes as `exp(-bx)/b`.
"""
function generate_trees(
    outdir; n=50, b=1, kwargs...
)
    n < 2 && error("Expected n>1, got n=$n")
    generate_trees(outdir, () -> genealogy(YuleCoalescent(n,b)); kwargs...)
end

function generate_trees(
    outdir, rand_tree::Function;
    M = 1,
    add_outgroup=true,
    outgroup_distance = :auto,
    outgroup_name = "outgroup",
)
    for r in 1:M
        out = outdir * "/$(r)/"
        mkpath(out)
        tree = rand_tree()
        # add outgroup
        if add_outgroup
            if outgroup_distance == :auto
                outgroup_distance = mean(l -> distance(tree.root, l), leaves(tree))
            end
            graft!(tree, TreeNode(;tau = outgroup_distance, label=outgroup_name), tree.root)
        end
        TreeTools.label!(tree, tree.root, "root")
        # write
        write(out * "/tree.nwk", tree)
    end
end

"""
    simulate_sequences(
        folder, model::ASR.EvolutionModel;
        tree_file, leaves_fasta, internals_fasta, prefix, kwargs...
    )

Use folder as base folder. Simulate sequences along `folder/tree_file` using `model`.
Store results as `folder/prefix/X_fasta` where `X` stands for leaves or internals.
Work is done by `ASRUtils.evolve`, additional `kwargs` are passed to it.
"""
function simulate_sequences(
    folder, model;
    tree_file = "tree.nwk",
    leaves_fasta = "alignment_leaves.fasta",
    internals_fasta = "alignment_internals.fasta",
    prefix = "",
    kwargs...
)
    tree = read_tree(folder * "/" * tree_file)
    return evolve(
        tree, model;
        leaves_fasta = folder * "/" * prefix *"/" * leaves_fasta,
        internals_fasta = folder * "/" * prefix *"/"* internals_fasta,
        kwargs...
    )
end
