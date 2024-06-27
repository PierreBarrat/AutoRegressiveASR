using ArDCA
using AutoRegressiveASR
using Chain
using CSV
using DataFrames
using DataFramesMeta
using DCATools
using DrWatson
using JSON3
using PlutoUI
using StatsBase
using StatsPlots

function sem(ystd, N)
    # error on mean calculated from sample standard deviation and number of samples
    # [-1.96, 196] has 95% of the mass for a normal distribution
    return 1.96 * ystd ./ sqrt.(N)
end

HAM = pluto_ingredients(
    scriptsdir("figures_and_results/analyze_results_and_write_df.jl")
)

function hamming_distance_plots(folder)

    # Read data

    data_all, _ = produce_or_load(
        Dict("folder" => folder);
        filename = x -> joinpath(x["folder"], "measures_asr.jld2"), suffix="",
    ) do config
        HAM.analyze_results_and_write(config["folder"])
    end
    data = data_all["asr"]

    simulation_parameters = JSON3.read(
        joinpath(folder, "simulation_parameters.json"), Dict
    )

    # Setup
    generative_model = @chain begin
        simulation_parameters["generative_model"]
        load_generative_model
    end

    model_consensus = let
        sample_file = projectdir(simulation_parameters["sample_equilibrium"])
        Seq = if isnothing(sample_file) || !isfile(sample_file)
            T = 100
            M = 1000
            @info "sampling generative model for consensus. Check eq. time (default $T)"
            S = DCATools.sample(generative_model, M; Twait = T)
        else
            read_msa(sample_file)
        end
        cons = DCATools.consensus(Seq) # a DCASample object
        DCATools.num_to_aa(cons[1], cons.mapping) # a string sequence
    end

    strategies = let
        st = collect(keys(data))

        lt(x,y) = if length(x) == length(y)
            x > y
        else
            length(x) > length(y)
        end
        sort(st; lt)
    end

    begin
        # smoothing alg
        w = 20
        outliers_right = 0.
        smoothing_alg = :hist
    end

    begin
        # plot style
        pal = palette(:default)
        strat_clr = Dict{Any,Any}(
            "iqtree" => pal[1], "autoregressive" => pal[2], "real" => pal[3]
        )
        for strat in strategies
            strat_clr[strat] = strat_clr[strat[1]]
        end
    end

    begin
        bayesian(strategies) = filter(x -> length(x)>1 && x[2]=="Bayes", strategies)
        ml(strategies) = filter(strategies) do x
            length(x) < 2 && return false
            x[2] == "ML" || x[2] == "ml"
        end
        real(strategies) = filter(==(("real",)), strategies)
        reconstruction(strategies) = filter(!=(("real",)), strategies)
        strat_label(strat) = joinpath(strat...)

        iqtree(strategies) = filter(x -> x[1]=="iqtree", strategies)
        ar(strategies) = filter(x -> x[1]=="autoregressive", strategies)

        function label_short(strat)
            length(strat) == 1 && return strat[1]
            strat[2] == "Bayes" ? "" : strat[1]
        end
        label_long(strat) = reduce((x,y) -> x*" - "*y, strat)


        function linestyle(strat)
            lw = 4
            return if length(strat) > 1 && strat[2] == "Bayes"
                (lw, :dash, strat_clr[strat[1]])
            else
                (lw, strat_clr[strat[1]])
            end
        end
        function barstyle(strat)
            (3, strat_clr[strat])
        end
    end

    # Plots

    hamming_real_ml_wgaps = let p = plot()
        for (i, strat) in enumerate(ml(strategies))
            x, y, ystd, N = ASRU.easy_smooth(
                data[strat], :node_depth, :hamming_to_real;
                w, alg=smoothing_alg, outliers_right
            )
            yerr = sem(ystd, N)
            plot!(
                x, y; ribbon = yerr, fillalpha=.2,
                label=label_short(strat), line=linestyle(strat)
            )
        end

        # Difference
        S1, S2 = (("iqtree", "ML"), ("autoregressive", "ML"))
        D1 = sort(data[S1], :node_depth)
        D2 = sort(data[S2], :node_depth)
        X = D1.node_depth
        Y = D1.hamming_to_real - D2.hamming_to_real # iqtree - AR

        x, y, ystd, N = ASRU.easy_smooth(
            X, Y; w, alg=smoothing_alg, outliers_right,
        )
        yerr = sem(ystd, N)
        plot!(
            x, y; ribbon = yerr, fillalpha=.2, label="improvement", color=:black
        )

        #
        plot!(
            xlabel = "Node depth",
            ylabel = "Hamming distance to real",
            title = "",
            frame = :box,
            legend = :topleft,
        )
        p
    end

    hamming_real_ml_nogaps = let p = plot()
        for (i, strat) in enumerate(ml(strategies))
            x, y, ystd, N = ASRU.easy_smooth(
                data[strat], :node_depth, :hamming_to_real_nogap;
                w, alg=smoothing_alg, outliers_right
            )
            yerr = sem(ystd, N)
            plot!(
                x, y; ribbon = yerr, fillalpha=.2,
                label=label_short(strat), line=linestyle(strat)
            )
        end

        # Difference
        S1, S2 = (("iqtree", "ML"), ("autoregressive", "ML"))
        D1 = sort(data[S1], :node_depth)
        D2 = sort(data[S2], :node_depth)
        X = D1.node_depth
        Y = D1.hamming_to_real_nogap - D2.hamming_to_real_nogap # iqtree - AR

        x, y, ystd, N = ASRU.easy_smooth(
            X, Y; w, alg=smoothing_alg, outliers_right,
        )
        yerr = sem(ystd, N)
        plot!(
            x, y; ribbon = yerr, fillalpha=.2, label="improvement", color=:black,
        )

        #
        plot!(
            xlabel = "Node depth",
            ylabel = "Hamming distance to real",
            title = "",
            frame = :box,
            legend = :topleft,
        )
        p
    end

    hamming_real_bayes_nogaps = let p = plot()
        for (i, strat) in enumerate(reconstruction(strategies))
            x, y, ystd, N = ASRU.easy_smooth(
                data[strat], :node_depth, :hamming_to_real_nogap;
                w, alg=smoothing_alg, outliers_right
            )
            if strat[2] == "Bayes"
                yerr = sem(ystd, N)
                plot!(
                    x, y; ribbon = yerr, fillalpha=.2,
                    label=strat[1], color = strat_clr[strat]
                )
            else
                plot!(
                    x, y; label="", color = strat_clr[strat], line = (3, :dash)
                )
            end
        end

        # Difference
        S1, S2 = (("iqtree", "Bayes"), ("autoregressive", "Bayes"))
        D1 = sort(data[S1], :node_depth)
        D2 = sort(data[S2], :node_depth)
        X = D1.node_depth
        Y = D1.hamming_to_real_nogap - D2.hamming_to_real_nogap # iqtree - AR

        x, y, ystd, N = ASRU.easy_smooth(
            X, Y; w, alg=smoothing_alg, outliers_right,
        )
        yerr = sem(ystd, N)
        plot!(
            x, y; ribbon = yerr, fillalpha=.2, label="improvement", color=:black
        )

        # Difference ML for ref
        S1, S2 = (("iqtree", "ML"), ("autoregressive", "ML"))
        D1 = sort(data[S1], :node_depth)
        D2 = sort(data[S2], :node_depth)
        X = D1.node_depth
        Y = D1.hamming_to_real_nogap - D2.hamming_to_real_nogap # iqtree - AR

        x, y, ystd, N = ASRU.easy_smooth(
            X, Y; w, alg=smoothing_alg, outliers_right,
        )
        plot!(
            x, y; label="", line = (:black, :dash, 3)
        )

        #
        plot!(
            xlabel = "Node depth",
            ylabel = "Hamming distance to real",
            title = "",
            frame = :box,
            legend = :topleft,
        )

        p
    end

    return (
        hamming_real_ml_wgaps = hamming_real_ml_wgaps,
        hamming_real_ml_nogaps = hamming_real_ml_nogaps,
        hamming_real_bayes_nogaps = hamming_real_bayes_nogaps,
    )
end
