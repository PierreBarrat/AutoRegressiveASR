function emd_hamming_to_ref(X::AbstractAlignment, ref::AbstractAlignment)
    step_x = div(length(X), 250)
    step_ref = div(length(ref), 250)
    Dx = pairwise_hamming(X, ref; step_left=step_x, step_right = step_ref)
    Dref = pairwise_hamming(ref; step = step_ref)

    normalize(x) = x / sum(x)
    h_edges = 0:.01:1
    Wx = fit(Histogram, vec(Dx), h_edges).weights |> normalize
    Wr = @chain [Dref[i,j] for i in 1:size(Dref,1) for j in (i+1):size(Dref,1)] begin
        fit(Histogram, _, h_edges).weights
        normalize
    end

    h_vals = (h_edges[2:end] + h_edges[1:end-1]) / 2
    # C = [(x - y)^2 for x in h_vals, y in h_vals]
    C = [abs(x-y) for x in h_vals, y in h_vals]

    return emd2(Wx, Wr, C, Tulip.Optimizer())
end

function emd_self_hamming(X::AbstractAlignment, Y::AbstractAlignment)
    step_x = div(length(X), 250)
    step_y = div(length(Y), 250)
    Dx = pairwise_hamming(X; step = step_x)
    Dy = pairwise_hamming(Y; step = step_y)

    normalize(x) = x / sum(x)
    h_edges = 0:.01:1
    Wx = @chain [Dx[i,j] for i in 1:size(Dx,1) for j in (i+1):size(Dx,1)] begin
        fit(Histogram, _, h_edges).weights
        normalize
    end
    Wy = @chain [Dy[i,j] for i in 1:size(Dy,1) for j in (i+1):size(Dy,1)] begin
        fit(Histogram, _, h_edges).weights
        normalize
    end

    h_vals = (h_edges[2:end] + h_edges[1:end-1]) / 2
    # C = [(x - y)^2 for x in h_vals, y in h_vals]
    C = [abs(x-y) for x in h_vals, y in h_vals]

    return emd2(Wx, Wy, C, Tulip.Optimizer())
end


function emd_cross_hamming(X::AbstractAlignment, Y::AbstractAlignment, Dref::AbstractVector)
    # compute cross-hamming between X and Y
    step_x = div(length(X), 250)
    step_y = div(length(Y), 250)
    Dxy = pairwise_hamming(X, Y; step = step_x, normalize = true) |> vec
    # fit histogram for Dxy and reference Dref
    h_edges = 0:.01:1
    h_vals, Wxy = fit_distance_histogram(h_edges, Dxy)
    _, Wref = fit_distance_histogram(h_edges, Dref)
    # compute EMD from the histograms
    C = [abs(x-y) for x in h_vals, y in h_vals]
    return emd2(Wxy, Wref, C, Tulip.Optimizer())
end

function emd_self_hamming(X::AbstractAlignment, Dref::AbstractVector)
    # Compute self hamming for X
    step_x = div(length(X), 250)
    Dx = pairwise_hamming(X; step = step_x, as_vec = true, normalize = true)
    # Fit histogram for Dx and reference
    h_edges = 0:.01:1
    h_vals, Wx = fit_distance_histogram(h_edges, Dx)
    _, Wref = fit_distance_histogram(h_edges, Dref)
    # Compute EMD from the histograms
    C = [abs(x-y) for x in h_vals, y in h_vals]
    return emd2(Wx, Wref, C, Tulip.Optimizer())
end

function fit_distance_histogram(h_edges, D::AbstractVector)
    normalize(x) = x / sum(x)
    W = fit(Histogram, D, h_edges).weights |> normalize
    h_vals = (h_edges[2:end] + h_edges[1:end-1]) / 2
    return h_vals, W
end
