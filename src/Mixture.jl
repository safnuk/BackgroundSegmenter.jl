import Base

const NOISE_SIGMA = 15.0
const INITIAL_VARIANCE = 4NOISE_SIGMA^2
# const MIN_VARIANCE = NOISE_SIGMA^2
const MIN_VARIANCE = 5.0^2
const BACKGROUND_RATIO = 0.7
const INITIAL_WEIGHT = .05
const VARIANCE_THRESHOLD = 2.5^2
const INITIALIZATION_WINDOW = 200
const MIN_LEARNING_RATE = 1.0 / INITIALIZATION_WINDOW

struct Gaussian
    μ :: Float64
    σ² :: Float64
    ω :: Float64
    sortkey :: Float64
    Gaussian(μ, σ², ω) = new(μ, σ², ω, ω / sqrt(σ²))
end

Gaussian(μ, σ²) = Gaussian(μ, σ², INITIAL_WEIGHT) 
Gaussian(μ) = Gaussian(μ, INITIAL_VARIANCE)
Base.isless(K::Gaussian, L::Gaussian) = K.sortkey < L.sortkey
neglog_likelihood(x, K::Gaussian) = (x - K.μ)^2 / K.σ²

Base.:*(K::Gaussian, x) = Gaussian(K.μ, K.σ², K.ω * x)
Base.:*(x, K::Gaussian) = K * x

mutable struct MixtureModel <: AbstractArray{Float64, 1}
    kernels :: Vector{Gaussian}
    threshold :: Float64
    background_threshold :: Float64
    number_frames :: Float64
    MixtureModel(kernels, threshold, background_threshold) = normalize!(new(kernels, threshold, background_threshold, 0.0))
end
function MixtureModel(kernels :: Vector{Gaussian})
    MixtureModel(kernels, VARIANCE_THRESHOLD, BACKGROUND_RATIO)
end
function MixtureModel(n::Int, μ=0, σ²=INITIAL_VARIANCE)
    MixtureModel([Gaussian(μ, σ²) for n in 1:n])
end
Base.size(m :: MixtureModel) = Base.size(m.kernels) 
Base.IndexStyle(::Type{MixtureModel}) = Base.IndexLinear()
Base.getindex(M::MixtureModel, i::Int) = Base.getindex(M.kernels, i)
Base.setindex!(M::MixtureModel, kernel::Gaussian, i::Int) = Base.setindex!(M.kernels, kernel, i)

function normalize!(M::MixtureModel; sort=true)
    total_weight = 0.0
    for kernel in M
        total_weight += kernel.ω
    end
    for (n, kernel) in enumerate(M)
        M[n] = kernel * (1 / total_weight)
    end
    if sort
        sort!(M.kernels, rev=true)
    end
    return M
end

function apply!(M::MixtureModel, x) :: UInt8
    k = 0
    matched = false
    M.number_frames += 1.0
    τ = M.threshold
    if M.number_frames > INITIALIZATION_WINDOW
        λ = MIN_LEARNING_RATE
    else
        λ = 1 / M.number_frames
    end
    for (k, kernel) in enumerate(M)
        if kernel.ω < eps() break end
        if neglog_likelihood(x, kernel) < τ
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
    if weight_sum < M.background_threshold
        return zero(UInt8)
    elseif x < average
        return zero(UInt8)
    else
        return one(UInt8)
    end
end

function bg_energy!(M::MixtureModel, x) :: Float64
    energy = Inf
    weight_sum = 0.0
    average = 0.0
    apply!(M, x)
    for kernel in M
        weight_sum += kernel.ω
        average += kernel.ω * kernel.μ
        energy = min(neglog_likelihood(x, kernel), energy)
        if weight_sum >= M.background_threshold
            break
        end
    end
    if x < average
        return min(energy, 0.9 * M.threshold)
    else
        return energy
    end
end

function update!(M::MixtureModel, matched_pos, x, λ)
    for (n, K) in enumerate(M)
        if n == matched_pos
            μ = K.μ + λ*(x - K.μ)
            σ² = max((1-λ)*K.σ² + λ*(x - K.μ)^2, MIN_VARIANCE)
            ω = (1-λ) * K.ω + λ
            M[n] =  Gaussian(μ, σ², ω )
        else
            M[n] = Gaussian(K.μ, K.σ², (1-λ)*K.ω )
        end
    end
    return move_kernel_into_order!(M, matched_pos)
end

function move_kernel_into_order!(M::MixtureModel, k)
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

function replace_kernel_at!(k::Int, M::MixtureModel, kernel::Gaussian)
    M[k] = kernel
    return move_kernel_into_order!(M, k)
end
