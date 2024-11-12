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

function hamming_distance_plots(folder; iqtree_strategy="base")
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
        lt(x,y) = if length(x) == length(y)
            x > y
        else
            length(x) > length(y)
        end
        function filter_iqtree(strat)
            if length(strat) != 2
                return true
            end

            name, _ = strat # iqtree-model / autoregressive
            return if occursin("iqtree", name)
                if iqtree_strategy == "base"
                    !occursin(r"C\d\d", name)
                else
                    occursin(Regex(iqtree_strategy), name) # select the wanted strategy
                end
            else
                true
            end
        end
        @info @chain data keys collect sort(; lt)
        st = @chain data keys collect sort(; lt) filter(filter_iqtree, _)
    end

    begin
        # smoothing alg
        w = 20
        outliers_right = 0.
        smoothing_alg = :hist
    end

    begin
        is_iqtree(strat) = occursin("iqtree", strat[1])
        is_autoregressive(strat) = strat[1] == "autoregressive"
        is_ml(strat) = length(strat) == 2 && (strat[2] == "ml" || strat[2] == "ML")
        is_bayes(strat) = length(strat) == 2 && (strat[2] == "Bayes" || strat[2] == "bayes")
    end

    begin
        pal = palette(:default)
        strat_clr = Dict()
        for strat in strategies
            if is_iqtree(strat)
                strat_clr[strat] = pal[1]
                strat_clr[strat[1]] = pal[1]
            elseif is_autoregressive(strat)
                strat_clr[strat] = pal[2]
                strat_clr[strat[1]] = pal[2]
            elseif strat[1] == "real"
                strat_clr[strat] = pal[3]
                strat_clr[strat[1]] = pal[3]
            end
        end
    end

    begin
        bayesian(strategies) = filter(is_bayes, strategies)
        ml(strategies) = filter(is_ml, strategies)
        real(strategies) = filter(==(("real",)), strategies)
        reconstruction(strategies) = filter(!=(("real",)), strategies)
        strat_label(strat) = joinpath(strat...)

        iqtree(strategies) = filter(is_iqtree , strategies)
        ar(strategies) = filter(is_autoregressive, strategies)


        function label_short(strat)
            length(strat) == 1 && return strat[1]
            strat[2] == "Bayes" ? "" : strat[1]
        end
        label_long(strat) = reduce((x,y) -> x*" - "*y, strat)


        function linestyle(strat)
            lw = 4
            return if is_bayes(strat)
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
        S_iqtree = @chain strategies filter(is_iqtree, _) filter(is_ml, _) first
        S_ar = @chain strategies filter(is_autoregressive, _) filter(is_ml, _) first
        D1 = sort(data[S_iqtree], :node_depth)
        D2 = sort(data[S_ar], :node_depth)
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
        S_iqtree = @chain strategies filter(is_iqtree, _) filter(is_ml, _) first
        S_ar = @chain strategies filter(is_autoregressive, _) filter(is_ml, _) first
        D1 = sort(data[S_iqtree], :node_depth)
        D2 = sort(data[S_ar], :node_depth)
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
        S_iqtree = @chain strategies filter(is_iqtree, _) filter(is_bayes, _) first
        S_ar = @chain strategies filter(is_autoregressive, _) filter(is_bayes, _) first
        D1 = sort(data[S_iqtree], :node_depth)
        D2 = sort(data[S_ar], :node_depth)
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
        S_iqtree = @chain strategies filter(is_iqtree, _) filter(is_ml, _) first
        S_ar = @chain strategies filter(is_autoregressive, _) filter(is_ml, _) first
        D1 = sort(data[S_iqtree], :node_depth)
        D2 = sort(data[S_ar], :node_depth)
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

    iqtree_strat_name = @chain strategies filter(is_iqtree, _) first first
    return (
        hamming_real_ml_wgaps = hamming_real_ml_wgaps,
        hamming_real_ml_nogaps = hamming_real_ml_nogaps,
        hamming_real_bayes_nogaps = hamming_real_bayes_nogaps,
    )
end
