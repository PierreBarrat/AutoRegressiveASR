using DrWatson
quickactivate(@__DIR__, "AutoRegressiveASR")

using AncestralSequenceReconstruction
using AutoRegressiveASR
using BioSequenceMappings
using Dates
using JLD2
using JSON3
using Random
using StatsBase
using TreeTools

const base_folder = datadir("Stiffler")
const ref_alignment_folder = datadir("Stiffler/aligned_data_ref/")
const R20_ref_alignment = joinpath(ref_alignment_folder, "PSE1_rnd20_aligned_PF13354_noinserts.fasta")
const wt_ref_alignment = joinpath(ref_alignment_folder, "PSE1_aligned_PF13354_noinserts.fasta")
const pfam_ref_alignment = joinpath(ref_alignment_folder, "PF13354_msa.fasta")

function make_tree(X::AbstractAlignment; time = .2)
    n = length(X)
    tree = star_tree(n, time)
    for (node, name) in zip(POTleaves(tree), X.names)
        label!(tree, node, name)
    end
    return tree
end

function select_data(M, params; kwargs...)
    return select_data(M, read_fasta(R20_ref_alignment), params; kwargs...)
end

function select_data(M, S, params; id = randstring(5))
    params = deepcopy(params)
    seqstyle = params["seq_style"]
    treestyle = params["tree_style"]

    outdir = joinpath(
        base_folder,
        "subalignments",
        "seqstyle_$(seqstyle)_treestyle_$(treestyle)$(params["suffix"])",
        "M$(M)",
        "id$(id)",
    )
    mkpath(outdir)
    outfile = joinpath(outdir, "R20.fasta")
    cp(wt_ref_alignment, joinpath(outdir, "wt.fasta"); force=true)
    mkpath(outdir)

    # Selecting sequences
    Su = if seqstyle == :uniform
        idx = randperm(sequence_number(S))[1:M]
        Su = subsample(S, idx)
    else
        error("Possible values for `style`: `(:uniform,)`")
    end
    write(outfile, Su)
    @info "Alignment written in $outfile"

    # Making tree
    tree = make_tree(Su)
    write(joinpath(outdir, "base_tree.nwk"), tree)

    # tag
    params["timestamp"] = now()
    @tag!(params)
    open(joinpath(outdir, "params.json"), "w") do f
        JSON3.pretty(f, JSON3.write(params))
    end

    return outdir, params
end

function reconstruct_all(
    folder, params;
    arnet = nothing, do_arnet=true, do_iqtree=true, do_cons=true
)
    tree_style = params["tree_style"]
    params = deepcopy(params)
    do_cons && (params["consensus"] = reconstruct_consensus(folder))
    do_iqtree && (params["iqtree"] = reconstruct_iqtree(
        folder; model = params["iqtree_model"], tree_style,
    ))
    do_arnet && (params["arnet"] = reconstruct_arnet(
        folder, arnet, params["arnet_file"];
        branch_opt = params["branch_opt"],
        tree_style,
        binarize = params["binarize"],
        with_code = params["with_code"],
    ))

    params["timestamp"] = now()
    @tag!(params)
    open(joinpath(folder, "params.json"), "w") do f
        JSON3.pretty(f, JSON3.write(params))
    end

    return params
end

function reconstruct_iqtree(
    folder;
    model="JTT+R6", scale_branch_length=true, root_with_init_tree=true, tree_style = :star,
)
    #= `folder` should contain
        - "wt.fasta"
        - "R20.fasta"
        - "base_tree.nwk"
    =#

    tree = joinpath(folder, "base_tree.nwk")
    aln = joinpath(folder, "R20.fasta")

    file_prefix = "IQTREE"
    mkpath(joinpath(folder, "iqtree"))
    tot_prefix = joinpath(folder, "iqtree", file_prefix)
    out_tree_file = "tree_iqtree.nwk"

    # reconstruction
    base_cmd = ["iqtree2 -asr -s $aln -pre $(tot_prefix) --keep-ident -redo"]
    if !isnothing(model) && !isempty(model)
        push!(base_cmd, "  -m $model ")
    end
    if tree_style == :star
        push!(base_cmd, " -te $tree ")
    end
    # what's below does not seem to work (iqtree problem?)
    # if scale_branch_length
    #     push!(base_cmd, " -blfix -blscale")
    # end

    cmd = string.(split(prod(base_cmd))) |> Cmd
    @info "Running iqtree with command $cmd"
    run(cmd)

    if root_with_init_tree
        tree_model = read_tree(tree)
        tree_iq = read_tree(tot_prefix*".treefile")
        tree_iq = ASRU.rearrange_from_model(tree_iq, tree_model)
        write(joinpath(dirname(tot_prefix), out_tree_file), tree_iq)
    else
        cp(
            joinpath(dirname(tot_prefix), "$(file_prefix).treefile"),
            joinpath(dirname(tot_prefix), out_tree_file);
            force=true
        )
    end
    cp(
        joinpath(dirname(tot_prefix), out_tree_file),
        joinpath(folder, out_tree_file)
    )

    return Dict("iqtree_model" => model)
end


function reconstruct_arnet(
    folder, arnet, arnet_file;
    branch_opt=:scale, tree_style = :star, binarize = true, with_code = false,
)
    @info "Reconstructing with arnet"
    outfolder = mkpath(joinpath(folder, "arnet_$(branch_opt)"))

    arnet = isnothing(arnet) ? JLD2.load(arnet_file)["arnet"] : arnet
    evo_arnet = ASR.AutoRegressiveModel(arnet; with_code)
    evo_profile = ASR.ProfileModel(arnet; with_code)

    tree_file = if tree_style == :star
        # binarized star for speed
        # !!! won't be compatible with branch_opt = :opt, but :star and :opt should not be called together
        if binarize
            tree = read_tree(joinpath(folder, "base_tree.nwk"))
            write(joinpath(folder, "base_tree_binarized.nwk"), binarize_star(tree))
            joinpath(folder, "base_tree_binarized.nwk")
        else
            @info "Not binarizing star tree"
            joinpath(folder, "base_tree.nwk")
        end
    else
        joinpath(folder, "tree_iqtree.nwk")
    end
    aln_file = joinpath(folder, "R20.fasta")
    cp(tree_file, joinpath(outfolder, "tree_start.nwk"))

    strategy_infer_bl = ASRMethod(; joint=false, optimize_branch_length_cycles=5, verbosity=2)
    strategy_ml = ASRMethod(; joint=false, ML=true, optimize_branch_length=false, verbosity=2)
    # strategy_bayes = ASRMethod(;
    #     joint=false, ML=false, repetitions = 10, optimize_branch_length=false, verbosity=2,
    # )

    # Recompute branch length
    tree_opt = joinpath(outfolder, "tree_opt.nwk")
    if branch_opt == :scale
        ASR.optimize_branch_scale(
            tree_file, aln_file, evo_profile, strategy_infer_bl;
            outnewick = tree_opt,
        )
    elseif branch_opt == :opt
        ASR.optimize_branch_length(
            tree_file, aln_file, evo_profile, strategy_infer_bl;
            outnewick = tree_opt,
        )
    end

    # ML inference
    infer_ancestral(
        tree_opt, aln_file, evo_arnet, strategy_ml;
        outfasta = joinpath(outfolder, "ML", "ml.fasta"),# outtable = joinpath(outfolder, "table.csv"), table_style = :verbose,
    )

    # # Bayesian inference
    # infer_ancestral(
    #     tree_opt, aln_file, evo_arnet, strategy_bayes;
    #     outfasta = joinpath(outfolder, "Bayes", "bayes.fasta"),
    #     alignment_per_node = true,
    # )

    return Dict(
        "asr_strategies" => Dict(
            "strat_bl" => strategy_infer_bl,
            "strat_ml" => strategy_ml,
            # "strat_bayes" => strategy_bayes,
        ),
        "arnet_file" => arnet_file,
        "tree_file" => tree_file,
    )
end

function binarize_star(tree)
    @assert length(nodes(tree)) - length(leaves(tree)) == 1 "Star tree should have only one internal node"
    T = branch_length(first(leaves(tree)))
    n = length(leaves(tree))
    btree = TreeTools.ladder_tree(n, missing)
    for node in nodes(btree)
        branch_length!(node, isleaf(node) ? T : 1e-6/n)
    end
    for (n1, n2) in zip(POTleaves(btree), POTleaves(tree))
        label!(btree, n1, label(n2))
    end
    return btree
end

function reconstruct_consensus(folder)
    @info "Reconstructing with consensus"
    aln = read_fasta(joinpath(folder, "R20.fasta"))
    mkpath(joinpath(folder, "consensus"))
    write(joinpath(folder, "consensus", "consensus.fasta"), consensus(aln))
    return nothing
end
