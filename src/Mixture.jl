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
Base.:/(K::Gaussian, x) = Gaussian(K.μ, K.σ², K.ω / x, K.sortkey / x)

struct MixtureParams
    threshold :: Float64
    background_threshold :: Float64
    learning_rate :: Float64
    number_frames :: Int
end

MixtureParams() = MixtureParams(VARIANCE_THRESHOLD, BACKGROUND_RATIO, 1.0, 0)

function increment_frames(p::MixtureParams)
    if p.number_frames > INITIALIZATION_WINDOW
        return MixtureParams(p.threshold, p.background_threshold,
                             MIN_LEARNING_RATE, p.number_frames+1)
    else
        return MixtureParams(p.threshold, p.background_threshold,
                             1.0 / (p.number_frames + 2), p.number_frames+1)
    end
end

mutable struct MixtureModel
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

function normalize!(M, idx)
    total_weight = 0.0
    n = size(M.kernels)[1]
    for i in 1:n
        total_weight += M.kernels[i, idx].ω
    end
    for i in 1:n
        M.kernels[i, idx] = M.kernels[i, idx] / total_weight
    end
    return M
end

function apply!(M::MixtureModel, input, output)
    for idx in CartesianRange(size(input))
        x = input[idx]
        output[idx] = _apply!(M, idx, x)
    end
    M.params = increment_frames(M.params)
    return output
end

function apply!(M::MixtureModel, input)
    output = similar(input)
    return apply!(M, input, output)
end

function _apply!(M::MixtureModel, idx, x)
    k = 0
    n = size(M.kernels)[1]
    matched = false
    τ = M.params.threshold
    λ = M.params.learning_rate
    for k in 1:n
        kernel = M.kernels[k, idx]
        if kernel.ω < eps() break end
        if neglog_likelihood(kernel, x) < τ
            k = update!(M, k, idx, x, λ)
            matched = true
            break
        end
    end
    if !matched
        kernel = Gaussian(x)
        k = replace_kernel_at!(k, idx, M, kernel)
    end
    weight_sum = 0.0
    average = 0.0
    for j in 1:k-1
        weight_sum += M.kernels[j, idx].ω
        average += M.kernels[j, idx].ω * M.kernels[j, idx].μ
    end
    if weight_sum < M.params.background_threshold
        return zero(x)
    elseif x < average
        return zero(x)
    else
        return one(x)
    end
end

function bg_energy!(M::MixtureModel, input)
    output = Array{Float64}(size(input))
    for idx in CartesianRange(size(input))
        x = input[idx]
        output[idx] = _bg_energy!(M, idx, x)
    end
    M.params = increment_frames(M.params)
    return output
end

function _bg_energy!(M, idx, x)
    energy = Inf
    weight_sum = 0.0
    average = 0.0
    _apply!(M, idx, x)
    n = size(M.kernels)[1]
    for i in 1:n
        kernel = M.kernels[i, idx]
        weight_sum += kernel.ω
        average += kernel.ω * kernel.μ
        energy = min(neglog_likelihood(kernel, x), energy)
        if weight_sum >= M.params.background_threshold
            break
        end
    end
    if x < average
        return min(energy, 0.9 * M.params.threshold)
    else
        return energy
    end
end

function update!(M, matched_pos, idx, x, λ)
    num_kernels = size(M.kernels)[1]
    for n in 1:num_kernels
        K = M.kernels[n, idx]
        if n == matched_pos
            μ = K.μ + λ*(x - K.μ)
            σ² = max((1-λ)*K.σ² + λ*(x - K.μ)^2, MIN_VARIANCE)
            ω = (1-λ) * K.ω + λ
            M.kernels[n, idx] =  Gaussian(μ, σ², ω)
        else
            M.kernels[n, idx] = (1-λ) * K
        end
    end
    return move_kernel_into_order!(M, matched_pos, idx)
end

function move_kernel_into_order!(M, k, idx)
    n = k
    for n in k-1:-1:1
        if M.kernels[n, idx].ω < M.kernels[n+1, idx].ω
            M.kernels[n, idx], M.kernels[n+1, idx] = M.kernels[n+1, idx], M.kernels[n, idx]
        else
            n = n+1
            break
        end
    end
    normalize!(M, idx)
    return n
end

function replace_kernel_at!(k, idx, M, kernel::Gaussian)
    M.kernels[k, idx] = kernel
    return move_kernel_into_order!(M, k, idx)
end
