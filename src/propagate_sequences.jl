#################################################################################
#################################### Potts ######################################
#################################################################################
"""
    propagate(s0::AbstractVector{Int}, steps::Int, g::DCAGraph[, M=1]; step_type=:metropolis)

Propagate Potts model `g` starting from `s0` and doing MCMC `steps`.
For Metropolis sampling, `steps` is the number of swaps. For Gibbs, it is the number of
sweeps.
"""
function propagate(s0::AbstractVector{Int}, steps::Int, g::DCAGraph; step_type=:metropolis)
    s = copy(s0)
    for _ in 1:steps
        _mcmc_step!(s, g; step_type)
    end
    return s
end
function propagate(
    s0::AbstractVector{Int}, steps::Int, g::DCAGraph, M::Int;
    step_type = :metropolis
)
    S = zeros(Int, M, length(s0))
    _s = copy(s0)
    for m in 1:M
        _s .= s0
        for _ in 1:steps
            _mcmc_step!(_s, g; step_type)
        end
        S[m, :] .= _s
    end
    return DCASample(S, g.q; mapping=g.mapping)
end

"""
    propagate(s0, times::AbstractVector{Int}, g::DCAGraph, M)

Propagate `g` starting from `s0`, sampling `M` sequences for each  at all times in `times`.
Return an array of `DCASample` corresponding to `times`.
"""
function propagate(s0::AbstractVector{Int}, times::AbstractVector{Int}, g::DCAGraph, M::Int)
    X = map(t -> zeros(Int, M, length(s0)), times)
    for m in 1:M
        S = _propagate(s0, times, g) # sample for each `t ∈ times`
        for t in eachindex(times)
            for i in 1:g.L
                X[t][m, i] = S[t, i]
            end
        end
    end
    return map(x -> DCASample(x, g.q; mapping=g.mapping), X)
end

"""
    _propagate(
        s0::AbstractVector{Int}, times::AbstractVector{Int}, g::DCAGraph;
        step_type = :metropolis
    )

single mcmc chain sampled at times in `times`.
"""
function _propagate(
    s0::AbstractVector{Int}, times::AbstractVector{Int}, g::DCAGraph;
    step_type = :metropolis
)
    S = zeros(Int, length(times), g.L)
    _s = copy(s0)
    (i1, t1) = (0, 0)
    while i1 < length(times)
        (i2, t2) = (i1+1, times[i1+1])
        Δt = t2 - t1
        for _ in 1:Δt
            _mcmc_step!(_s, g; step_type)
        end
        S[i2, :] .= _s

        (i1, t1) = (i2, t2)
    end
    return S
end

function _mcmc_step!(s, g::DCAGraph; step_type = :metropolis)
    if step_type == :metropolis
        DCATools.metropolis_step!(s, g)
    elseif step_type == :gibbs
        DCATools.gibbs_step!(s, g)
    else
        throw(ErrorException("Unrecognized step type: pick from `:metropolis` or `:gibbs`"))
    end
end


#################################################################################
############################### AutoRegressive ##################################
#################################################################################

"""
    propagate(s0, t::Number, g::ArNet, M::Int)
    propagate(s0, tvals::Array, g::ArNet, M::Int)

Sample `M` sequences at time `t` from `s0`.
"""
function propagate(s0, t::Number, g::ArNet, M::Int)
    q = length(g.p0)
    L = length(g.H) # it's L-1
    backorder = sortperm(g.idxperm)
    s0_perm = s0[g.idxperm]

    res = Matrix{Int}(undef, L+1, M) # sample
    totH = Vector{Float64}(undef, q) # local field
    for m in 1:M
        res[1, m] = rand() < exp(-t) ? s0_perm[1] : Distributions.wsample(1:q, g.p0)
        for site in 1:L
            if rand() < exp(-t)
                res[site+1, m] = s0_perm[site+1]
            else
                copy!(totH, g.H[site])
                for i in 1:site, a in 1:q
                    totH[a] += g.J[site][a, res[i, m], i]
                end
                ArDCA.softmax!(totH)
                res[site+1, m] = Distributions.wsample(1:q, totH)
            end
        end
    end
    ArDCA.permuterow!(res, backorder)
    return DCASample(res', q)
end

propagate(s0, tvals::AbstractVector, g::ArNet, M) = map(t -> propagate(s0, t, g, M), tvals)
propagate(s0, t::Number, g::ArNet) = propagate(s0, t, g, 1)[1]

#################################################################################
############################# DCATools.ProfileModel #############################
#################################################################################

"""
    propagate(s0, t::Number, g::DCATools.ProfileModel, M::Int)
    propagate(s0, tvals::Array, g::DCATools.ProfileModel, M::Int)

Sample `M` sequences at time `t` from `s0`. `DCATools.ProfileModel` is from `DCATools`.
"""
function propagate(s0, t::Number, g::DCATools.ProfileModel, M::Int)
    q, L = size(g)
    ν = exp(-t)
    res = Matrix{Int}(undef, L, M)
    for m in 1:M, i in 1:L
        res[i,m] = rand() < ν ? s0[i] : Distributions.wsample(1:q, g[i, :])
    end
    return DCASample(res', q)
end
function propagate(s0, tvals::AbstractVector, g::DCATools.ProfileModel, M)
    return map(t -> propagate(s0, t, g, M), tvals)
end
propagate(s0, t::Number, g::DCATools.ProfileModel) = propagate(s0, t, g, 1)[1]
