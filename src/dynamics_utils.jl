function useful_callbacks()
    h_to_i(S, model, s0, t; kwargs...) = map(s -> DCATools.hamming(s, s0), S)
    en(S, model, s0, t; ref_model=model) = map(s -> energy(ref_model, s), S)
    return (
        energy = en,
        hamming_to_init = h_to_i,
        pw_hamming = (S, model, s0, t; kwargs...) -> DCATools.pw_hamming_distance(S; normalize=false),
    )
end

function average_at_t(
    s0::AbstractVector{<:Int}, functions, model, tvals, M=100; kwargs...
)
    time_samples = propagate(s0, tvals, model, M)
    return map(zip(tvals, time_samples)) do (t, S)
        names = (:t, Iterators.filter(!=(:t), keys(functions))...)
        values = (
            t,
            (@chain functions begin
                pairs
                Iterators.filter(p -> p[1] != :t, _)
                map(p -> p[2](S, model, s0, t; kwargs...), _)
                map(mean, _)
            end)...
        )
        NamedTuple{names}(values)
    end
end

"""
    average_energy_trajectory(S0::DCASample, g, potts, tvals; M=100)

For each `s` in `S0`, return <H(t)> for `t` in `tvals` and where the average runs over
a sample of `M` sequences. Uses `g` as a model for propagation and `potts` to
compute the energy.
"""
function average_energy_trajectory(
    S0::DCASample, arnet::ArNet, potts::DCAGraph, tvals;
    M = 100, scale = true,
)
    E = map(S0) do s0
        average_energy_trajectory(s0, arnet, potts, tvals; M)
    end
    Emean = mean(s -> energy(potts, s), DCATools.sample(potts, 5*M, Twait = 5*potts.L))
    return if scale
        map(x -> x .- Emean, E)
    else
        E
    end
end

function average_energy_trajectory(
    s0::AbstractVector{<:Int}, arnet::ArNet, potts, tvals; M = 100
)
    return map(tvals) do t
        S = propagate(s0, t, arnet, M)
        mean(s -> energy(potts, s), S)
    end
end

function ground_state(arnet::ArNet)
    q = length(arnet.p0)
    L = length(arnet.H) + 1
    holder = Vector{Float64}(undef, q)

    x = Vector{Int64}(undef, L)
    for i in 1:L
        cond_prob = _conditional_probability(x, i, arnet, holder)
        x[i] = argmax(cond_prob)
    end

    return x[sortperm(arnet.idxperm)]
    # return x
end
