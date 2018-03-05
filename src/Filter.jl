struct Filter
    w::Array{Float64, 3}
    n::Int
    threshold::Float64
    Filter(n, t) = new(weights(n), n, t)
end
Filter(n) = Filter(n, 0.5)

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
                              temp::AbstractArray{T, 2},
                              w::AbstractArray{T, 2}, n=1) where {T}
    @assert size(out) == size(w) == size(temp)
    (a, b) = size(w)
    n = 3
    for j in 1:b, i in 1:a
        if w[i, j] == zero(T) 
            continue 
        end
        h_start, h_end = calc_bounds(i, a, n)
        w_start, w_end = calc_bounds(j, b, n)
        for q in w_start:w_end, p in h_start:h_end
            temp[i+p, j+q] = one(T)
        end
    end
    for j in 1:b, i in 1:a
        if temp[i, j] == zero(T) 
            continue 
        end
        out[i, j] = one(T)
        h_start, h_end = calc_bounds(i, a, n)
        w_start, w_end = calc_bounds(j, b, n)
        for q in w_start:w_end, p in h_start:h_end
            if temp[i+p, j+q] == zero(T)
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
