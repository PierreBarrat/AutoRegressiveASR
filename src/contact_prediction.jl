function ppv(true_distances::DataFrame, scores::DataFrame; Δi = 4, contact_threshold=8)
    @assert names(true_distances) == ["i", "j", "distance_A"] "wrong df format $(names(true_distances))"
    @assert names(scores) == ["i", "j", "score"] "wrong dataframe format $(names(scores))"
    dat = @chain begin
        innerjoin(true_distances, scores; on = [:i, :j])
        @subset abs.(:i - :j) .> Δi
        @transform @byrow :contact = (:distance_A < contact_threshold)
        @orderby(-:score, rev=true)
    end
    return cumsum(dat.contact) ./ collect(1:size(dat,1))
end

function ppv(true_distances::AbstractString, scores::AbstractString; kwargs...)
    return ppv(
        DataFrame(CSV.File(true_distances)),
        DataFrame(CSV.File(scores));
        kwargs...
    )
end
