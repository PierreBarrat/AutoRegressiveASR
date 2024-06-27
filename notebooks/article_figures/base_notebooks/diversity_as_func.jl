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

DIV = pluto_ingredients(scriptsdir("figures_and_results/diversity_functions.jl"))

function sem(ystd, N)
    # error on mean calculated from sample standard deviation and number of samples
    # [-1.96, 196] has 95% of the mass for a normal distribution
    return 1.96 * ystd ./ sqrt.(N)
end

function diversity_plots(folder)


    path_folder, data_folder, prefix = let
        p, d = dirname(folder), basename(folder)
        prefix = split(d, "_")[1]
        p, d, prefix
    end


    data, filename = produce_or_load(
       Dict("basefolder" => folder, "out" => "diversity_data.jld2");
       filename = x -> joinpath(x["basefolder"], x["out"]),
       suffix = "",
    ) do config
       data = DIV.diversity_data(config["basefolder"], config["out"])
    end

    strategies = ["iqtree", "ardca"]

    begin
        # smoothing width
        w = 20
        outliers_right = 0.
        smoothing_alg = :hist
        pal = palette(:default)
        strat_color = Dict(strat => pal[i] for (i, strat) in enumerate(strategies))
    end

    strat_label = Dict(
        "iqtree" => "iqtree",
        "ardca" => "autoregressive",
    )

    linestyle = let
        lw = 4
        Dict(strat => (lw, strat_color[strat]) for strat in strategies)
    end

    # Figures
    selfhamming = let p = plot()
        for strat in strategies
            x, y, ystd, N = ASRU.easy_smooth(
                data[strat], :depth, :av_self_hamming;
                w, alg=smoothing_alg, outliers_right,
            )
            yerr = sem(ystd, N)
            plot!(
                x, y, ribbon = yerr;
                fillalpha = .2, label=strat_label[strat], line=linestyle[strat]
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

    entropy = let p = plot()
        for strat in strategies
            x, y = ASRU.easy_smooth(
                data[strat], :depth, :entropy; w, alg=smoothing_alg, outliers_right,
            )
            plot!(x, y, label=strat_label[strat], line=linestyle[strat])
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

    return (selfhamming = selfhamming, entropy = entropy)
end
