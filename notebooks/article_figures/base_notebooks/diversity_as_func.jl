quickactivate(@__DIR__, "AutoRegressiveASR")

using ArDCA
using AutoRegressiveASR
using CSV
using DataFrames
using DataFramesMeta
using DCATools
using DrWatson
using JSON3
using PlutoUI
using StatsBase
using StatsPlots

includet(scriptsdir("figures_and_results/diversity_functions.jl"))

function _sem(ystd, N)
    # error on mean calculated from sample standard deviation and number of samples
    # [-1.96, 196] has 95% of the mass for a normal distribution
    return 1.96 * ystd ./ sqrt.(N)
end

function diversity_plots(folder; iqtree_strategy="base")
    @info "Diversity plots for $folder"
    dat, filename = produce_or_load(
       Dict("basefolder" => folder, "out" => "diversity_data.jld2");
       filename = x -> joinpath(x["basefolder"], x["out"]),
       suffix = "",
    ) do config
       dat = diversity_data(config["basefolder"], config["out"])
    end
    dat = dat["diversity"]

    strategies = let
        lt(x,y) = if length(x) == length(y)
            x > y
        else
            length(x) > length(y)
        end
        function filter_iqtree(strat)
            return if occursin("iqtree", strat)
                if iqtree_strategy == "base"
                    !occursin(r"C\d\d", strat)
                else
                    occursin(Regex(iqtree_strategy), strat) # select the wanted strategy
                end
            else
                true
            end
        end
        st = @chain dat keys collect sort(; lt) filter(filter_iqtree, _)
    end

    begin
        # smoothing width
        w = 20
        outliers_right = 0.
        smoothing_alg = :hist
    end
    begin
        is_iqtree(strat) = occursin("iqtree", strat)
        is_autoregressive(strat) = strat == "autoregressive"
    end
    begin
        pal = palette(:default)
        strat_clr = Dict()
        for strat in strategies
            if is_iqtree(strat)
                strat_clr[strat] = pal[1]
            elseif is_autoregressive(strat)
                strat_clr[strat] = pal[2]
            end
        end
    end

    linestyle = let
        lw = 4
        Dict(strat => (lw, strat_clr[strat]) for strat in strategies)
    end

    # Figures
    selfhamming = let p = plot()
        for strat in strategies
            x, y, ystd, N = ASRU.easy_smooth(
                dat[strat], :depth, :av_self_hamming;
                w, alg=smoothing_alg, outliers_right,
            )
            yerr = _sem(ystd, N)
            plot!(
                x, y, ribbon = yerr;
                fillalpha = .2, label=strat, line=linestyle[strat], color=strat_clr[strat]
            )
        end
        plot!(
            xlabel = "Node depth",
            ylabel = "Self-Hamming distance",
            # xlim = (-0.025, 2.025),
            title = "",
            frame = :box,
            legend = :topleft,
        )
        p
    end

    plt_entropy = let p = plot()
        for strat in strategies
            x, y = ASRU.easy_smooth(
                dat[strat], :depth, :entropy; w, alg=smoothing_alg, outliers_right,
            )
            plot!(x, y, label=strat, line=linestyle[strat], color=strat_clr[strat])
        end
        plot!(
            xlabel = "Node depth",
            ylabel = "Self-Hamming distance",
            title = "",
            frame = :box,
            legend = :topleft,
        )
        p
    end

    return (selfhamming = selfhamming, entropy = plt_entropy)
end
