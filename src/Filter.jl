using Clustering

struct Filter
    w::Array{Float64, 3}
    n::Int
    threshold::Float64
    Filter(n, t) = new(weights(n), n, t)
end
Filter(n) = Filter(n, 0.5)

function cluster(w::AbstractArray{T, 3}, radius=3; min_neighbors=20, min_cluster_size=30) where {T}
    (n, m, t) = size(w)
    fg_clusters = zeros(w)
    fg_points = hcat(map(x -> collect(convert.(Float64, ind2sub(w, x))), find(w))...)
    C = dbscan(fg_points, radius; min_neighbors=min_neighbors, min_cluster_size=min_cluster_size)
    for x in C
        for y in x.core_indices 
            index = convert.(Int, fg_points[:, y])
            fg_clusters[index...] = 1
        end
        # for y in x.boundary_indices 
        #     index = convert.(Int, fg_points[:, y])
        #     video[index...] = colorant"red"
        # end
    end
    return fg_clusters
end

function weights(n)
    w = zeros(Float64, 2n+1, 2n+1, 2n+1)
    for i in 1:2n+1, j in 1:2n+1, k in 1:2n+1
        if n * sqrt(3) > sqrt((i-n-1)^2 + (j-n-1)^2 + (k-n-1)^2)
            w[i, j, k] = 1.0
        end
        # w[i, j, k] = n * sqrt(3) - sqrt((i-n-1)^2 + (j-n-1)^2 + (k-n-1)^2)
    end
    w .= w ./ sum(w)
    return w
end

function morphological_close(w::AbstractArray{T, 2}, n=1) where {T}
    out = zeros(w)
    temp = zeros(w)
    morphological_close!(out, temp, w, n)
end

function morphological_close!(out::AbstractArray{T, 2},
                              temp::AbstractArray{S, 2},
                              w::AbstractArray{R, 2}, n=1) where {R, S, T}
    @assert size(out) == size(w) == size(temp)
    (a, b) = size(w)
    n = 3
    m = 6
    for j in 1:b, i in 1:a
        if w[i, j] == zero(R) 
            continue 
        end
        h_start, h_end = calc_bounds(i, a, n)
        w_start, w_end = calc_bounds(j, b, m)
        for q in w_start:w_end, p in h_start:h_end
            temp[i+p, j+q] = one(S)
        end
    end
    for j in 1:b, i in 1:a
        if temp[i, j] == zero(S) 
            continue 
        end
        out[i, j] = one(T)
        h_start, h_end = calc_bounds(i, a, n)
        w_start, w_end = calc_bounds(j, b, m)
        for q in w_start:w_end, p in h_start:h_end
            if temp[i+p, j+q] == zero(S)
                out[i, j] = zero(T)
                break
            end
        end
    end
    return out
end


function apply(f::Filter, w::Array{T, 3}) where {T}
    result = zeros(Float64, size(w))
    (a, b, c) = size(w)
    n = f.n
    for k in 1:c, j in 1:b, i in 1:a
        if w[i, j, k] == zero(T) 
            continue 
        end
        h_start, h_end = calc_bounds(i, a, n)
        w_start, w_end = calc_bounds(j, b, n)
        d_start, d_end = calc_bounds(k, c, n)
        for r in d_start:d_end, q in w_start:w_end, p in h_start:h_end
            result[i+p, j+q, k+r] += f.w[n+1+p, n+1+q, n+1+r]
        end
    end
    for (idx, x) in enumerate(result)
        if x < f.threshold
            result[idx] = 0.0
        else
            result[idx] = 1.0
        end
    end
    return result
end

function scale(w::Array{UInt8, 3}, s::Array{Float64, 3})
    @assert size(w) == size(s)
    result = zeros(w)
    for (idx, x) in enumerate(w)
        result[idx] = convert(UInt8, round(float(x) * s[idx]))
    end
    return result
end

function calc_bounds(idx, upper_limit, filter_size)
    return (max(-filter_size, 1 - idx), min(filter_size, upper_limit - idx))
end
