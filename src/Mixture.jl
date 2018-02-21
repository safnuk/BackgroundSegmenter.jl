import Base
using StaticArrays: SVector

const NOISE_SIGMA = 15.0
const INITIAL_VARIANCE = 4NOISE_SIGMA^2
# const MIN_VARIANCE = NOISE_SIGMA^2
const MIN_VARIANCE = 5.0^2
const BACKGROUND_RATIO = 0.7
const INITIAL_WEIGHT = .05
const INITIAL_SORTKEY = INITIAL_WEIGHT / sqrt(INITIAL_VARIANCE)
const VARIANCE_THRESHOLD = 2.5^2
const INITIALIZATION_WINDOW = 200
const MIN_LEARNING_RATE = 1.0 / INITIALIZATION_WINDOW

struct Gaussian
    μ :: Float64
    σ² :: Float64
    ω :: Float64
    sortkey :: Float64
end
Gaussian(μ, σ², ω) = Gaussian(μ, σ², ω, ω / sqrt(σ²))
Gaussian(μ, σ²) = Gaussian(μ, σ², INITIAL_WEIGHT)
Gaussian(μ) = Gaussian(μ, INITIAL_VARIANCE, INITIAL_WEIGHT, INITIAL_SORTKEY)
Base.isless(K::Gaussian, L::Gaussian) = K.sortkey < L.sortkey
neglog_likelihood(K::Gaussian, x) = (x - K.μ)^2 / K.σ²

Base.:*(K::Gaussian, x) = Gaussian(K.μ, K.σ², K.ω * x, K.sortkey * x)
Base.:*(x, K::Gaussian) = K * x

mutable struct MixtureParams
    threshold :: Float64
    background_threshold :: Float64
    learning_rate :: Float64
    number_frames :: Int
end

MixtureParams() = MixtureParams(VARIANCE_THRESHOLD, BACKGROUND_RATIO, 1.0, 0)

function increment_frames!(p::MixtureParams)
    p.number_frames += 1
    if p.number_frames > INITIALIZATION_WINDOW
        p.learning_rate = MIN_LEARNING_RATE
    else
        p.learning_rate = 1 / (p.number_frames + 1)
    end
    return p
end

struct MixtureModel
    kernels::Array{Gaussian, 3}
    n::Int
    params::MixtureParams
end

function MixtureModel(number_gaussians::Int, rows::Int, cols::Int, μ=0.0, σ²=INITIAL_VARIANCE)
    kernels = Array{Gaussian}(number_gaussians, rows, cols)
    weight = 1.0 / number_gaussians
    kernels .= Gaussian(μ, σ², weight)
    MixtureModel(kernels, number_gaussians, MixtureParams())
end
MixtureModel(n::Int, m::Int) = MixtureModel(5, n, m)

#Base.size(m :: MixtureModel) = Base.size(m.kernels)
#Base.IndexStyle(::Type{MixtureModel}) = Base.IndexLinear()
#Base.getindex(M::MixtureModel, i::Int) = Base.getindex(M.kernels, i)
#Base.setindex!(M::MixtureModel, kernel::Gaussian, i::Int) = Base.setindex!(M.kernels, kernel, i)

function normalize!(model; sort=true)
    total_weight = 0.0
    for kernel in model
        total_weight += kernel.ω
    end
    for (i, kernel) in enumerate(model)
        model[i] = kernel * (1 / total_weight)
    end
    if sort
        sort!(model, rev=true)
    end
    return model
end

function apply!(M::MixtureModel, input, output)
    (rows, cols) = size(input)
    for m in 1:cols, n in 1:rows
        output[n, m] = _apply!((@view (M.kernels[:, n, m])), M.params, input[n, m])
    end
    increment_frames!(M.params)
    return output
end

function apply!(M::MixtureModel, input)
    output = similar(input)
    return apply!(M, input, output)
end

function _apply!(M, params::MixtureParams, x)
    k = 0
    matched = false
    τ = params.threshold
    λ = params.learning_rate
    for (k, kernel) in enumerate(M)
        if kernel.ω < eps() break end
        if neglog_likelihood(kernel, x) < τ
            k = update!(M, k, x, λ)
            matched = true
            break
        end
    end
    if !matched
        kernel = Gaussian(x)
        k = replace_kernel_at!(k, M, kernel)
    end
    weight_sum = 0.0
    average = 0.0
    for n in 1:k-1
        weight_sum += M[n].ω
        average += M[n].ω * M[n].μ
    end
    if weight_sum < params.background_threshold
        return zero(x)
    elseif x < average
        return zero(x)
    else
        return one(x)
    end
end

function bg_energy!(M::MixtureModel, x)
    result = Array{Float64}(size(x))
    for (n, z) in enumerate(x)
        result[n] = _bg_energy!((@view M.kernels[:, n]), M.params, z)
    end
    increment_frames!(M.params)
    return result
end

function _bg_energy!(M, params::MixtureParams, x)
    energy = Inf
    weight_sum = 0.0
    average = 0.0
    _apply!(M, params, x)
    for kernel in M
        weight_sum += kernel.ω
        average += kernel.ω * kernel.μ
        energy = min(neglog_likelihood(kernel, x), energy)
        if weight_sum >= params.background_threshold
            break
        end
    end
    if x < average
        return min(energy, 0.9 * params.threshold)
    else
        return energy
    end
end

function update!(M, matched_pos, x, λ)
    for (n, K) in enumerate(M)
        if n == matched_pos
            μ = K.μ + λ*(x - K.μ)
            σ² = max((1-λ)*K.σ² + λ*(x - K.μ)^2, MIN_VARIANCE)
            ω = (1-λ) * K.ω + λ
            M[n] =  Gaussian(μ, σ², ω)
        else
            M[n] = (1-λ) * K
        end
    end
    return move_kernel_into_order!(M, matched_pos)
end

function move_kernel_into_order!(M, k)
    n = k
    for n in k-1:-1:1
        if M[n].ω < M[n+1].ω
            M[n], M[n+1] = M[n+1], M[n]
        else
            n = n+1
            break
        end
    end
    normalize!(M, sort=false)
    return n
end

function replace_kernel_at!(k::Int, M, kernel::Gaussian)
    M[k] = kernel
    return move_kernel_into_order!(M, k)
end
