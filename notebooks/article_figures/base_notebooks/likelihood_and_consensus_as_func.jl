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

function likelihood_and_hamming_consensus(folder)
    # Setup
    data_all, _ = produce_or_load(
        Dict("folder" => folder);
        filename = x -> joinpath(x["folder"], "measures_asr.jld2"), suffix="",
    ) do config
        HAM.analyze_results_and_write(config["folder"])
    end;

    data = data_all["asr"];

    simulation_parameters = JSON3.read(
        joinpath(folder, "simulation_parameters.json"), Dict
    )

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
    likelihood_v_depth = let p = plot()
        for strat in vcat(ml(strategies), bayesian(strategies), real(strategies))
            x, y = ASRU.easy_smooth(
                data[strat], :node_depth, :loglikelihood;
                w, alg=smoothing_alg, outliers_right,
            )
            plot!(x, y, label=label_short(strat), line=linestyle(strat))
        end

        # hline!([mean(data_all["likelihood_eq"])], line = (2, :black, :dash), label="")

        lk_lim = let
            r, l = (1.15, 4)
            xmin, xmax = extrema(data_all["likelihood_eq"])
            mid = (xmax + xmin)/2
            mid - (mid - xmin)/l, mid + (xmax - mid)/r
        end
        density_lk_eq = let
            plt = density(data_all["likelihood_eq"]; xlim=lk_lim)
            x, y = plt[1][1][:x], plt[1][1][:y]
            plot(
                y, x;
                ylim=lk_lim, label="", xticks = [], yticks = [], line=(3, :black),
                frame=:no, xaxis=false, fill=true, fillcolor = :black, fillalpha=.1
            )
        end

        plot!(p;
            xlabel = "Node depth",
            ylabel = "Likelihood",
            # xlim = (-0.025, 2.025),
            frame = :box,
            legend = :topleft,
            ylim = lk_lim,
        )

        p
        density_lk_eq

        plot(
            p, density_lk_eq;
            layout = @layout [a{0.95w} _ b{0.2w}]
        )
    end

    hamming_to_consensus = let p = plot()
        for strat in strategies
            x, y = ASRU.easy_smooth(
                data[strat], :node_depth, :hamming_to_aln_consensus;
                w, alg=smoothing_alg, outliers_right,
            )
            plot!(x, y, label=label_short(strat), line=linestyle(strat))
        end

        #

        plot!(
            xlabel = "Node depth",
            # xlim = (-0.025, 2.025),
            ylabel = "Hamming distance to consensus",
            frame = :box,
            # legend = :bottomleft,
        )

        p
    end

    return (
        likelihood_v_depth = likelihood_v_depth,
        hamming_to_consensus = hamming_to_consensus
    )

end
