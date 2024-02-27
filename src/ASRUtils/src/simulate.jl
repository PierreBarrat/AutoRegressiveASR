# """
#     generate_trees(outdir; n=50, b=1, M=10)

# For `r` in `1:M` generaet a hierarchy of folders:
# ```
# outdir/r/tree.nwk
# ```
# with `tree` sampled from the Yule coalescent with parameter `b` and `n` lineages.
# The average height of trees obtained like this is `log(n)/b`.
# The number of internal nodes distant of `x` from the leaves goes as `exp(-bx)/b`.
# """
# function generate_trees(
#     outdir; n=50, b=1, kwargs...
# )
#     n < 2 && error("Expected n>1, got n=$n")
#     generate_trees(outdir, () -> genealogy(YuleCoalescent(n,b)); kwargs...)
# end

# function generate_trees(
#     outdir, rand_tree::Function;
#     M = 1,
#     add_outgroup=true,
#     outgroup_distance = :auto,
#     outgroup_name = "outgroup",
# )
#     for r in 1:M
#         out = outdir * "/$(r)/"
#         mkpath(out)
#         tree = rand_tree()
#         # add outgroup
#         if add_outgroup
#             if outgroup_distance == :auto
#                 outgroup_distance = mean(l -> distance(tree.root, l), leaves(tree))
#             end
#             graft!(tree, TreeNode(;tau = outgroup_distance, label=outgroup_name), tree.root)
#         end
#         label!(tree, tree.root, "root")
#         # write
#         write(out * "/tree.nwk", tree)
#     end
# end

# """
#     simulate_sequences(
#         folder, model::ASR.EvolutionModel;
#         tree_file, leaves_fasta, internals_fasta, prefix, kwargs...
#     )

# Use folder as base folder. Simulate sequences along `folder/tree_file` using `model`.
# Store results as `folder/prefix/X_fasta` where `X` stands for leaves or internals.
# Work is done by `ASRUtils.evolve`, additional `kwargs` are passed to it.
# """
# function simulate_sequences(
#     folder, model;
#     tree_file = "tree.nwk",
#     leaves_fasta = "alignment_leaves.fasta",
#     internals_fasta = "alignment_internals.fasta",
#     prefix = "",
#     kwargs...
# )
#     tree = read_tree(folder * "/" * tree_file)
#     return evolve(
#         tree, model;
#         leaves_fasta = folder * "/" * prefix *"/" * leaves_fasta,
#         internals_fasta = folder * "/" * prefix *"/"* internals_fasta,
#         kwargs...
#     )
# end

# function simulate_sequences(
#     folder, arnet::ArDCA.ArNet;
#     tree_file = "tree.nwk",
#     leaves_fasta = "alignment_leaves.fasta",
#     internals_fasta = "alignment_internals.fasta",
#     prefix = "",
#     kwargs...
# )
# # this uses DCA_dynamics_tools.evolve_tree
#     tree = read_tree(folder * "/" * tree_file)
#     return evolve(
#         tree, arnet;
#         leaves_fasta = folder * "/" * prefix *"/" * leaves_fasta,
#         internals_fasta = folder * "/" * prefix *"/"* internals_fasta,
#         kwargs...
#     )
# end



"""
    infer_branch_length(
        folder, model::ASR.ProfileModel;
        tree_file, alignment_file, prefix, out_tree_file,
    )

Recompute branches of `folder/tree_file` using `model`
and sequences `folder/alignment_file`.
Output tree is written to `folder/prefix/out_tree_file`.
"""
function infer_branch_length(
    folder, model::ASR.ProfileModel;
    tree_file = "tree.nwk",
    alignment_file = "alignment_leaves.fasta",
    prefix = "",
    out_tree_file = "tree_recomputed_branches.nwk",
    alphabet = :aa,
)
    mkpath(joinpath(folder, prefix))
    strategy = ASR.ASRMethod(; joint=false, alphabet)
    return ASR.optimize_branch_length(
        folder * "/" * tree_file,
        folder * "/" * alignment_file,
        model,
        strategy;
        outnewick = joinpath(folder, prefix, out_tree_file)
    )
end



function infer_branch_length_iqtree(
    folder;
    tree_file = "tree.nwk",
    alignment_file = "alignment_leaves.fasta",
    prefix = "",
    iqtree_prefix = "IQTREE",
    outgroup_name = "outgroup",
    out_tree_file = "tree_recomputed_branches.nwk",
)

    tree = joinpath(folder, tree_file)
    aln = joinpath(folder, alignment_file)
    mkpath(joinpath(folder, prefix))
    tot_prefix = joinpath(folder, prefix, iqtree_prefix)
    run(`iqtree2 -s $aln -te $tree -pre $(tot_prefix) -o $(outgroup_name) --keep-ident -redo`)
    # copying IQTREE.treefile (or )
    cp(
        joinpath(folder, prefix, "$(iqtree_prefix).treefile"),
        joinpath(folder, prefix, out_tree_file)
    )
end
