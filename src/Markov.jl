struct MarkovModel
    background::Array{MixtureModel, 2}
    params::MixtureParams
    cut::MinCut
    function MarkovModel(G, kernels=5, smoothing=2.0)
        (n, m) = size(G)
        cut = MinCut(zeros(Float64, n, m),
                     VARIANCE_THRESHOLD,
                     smoothing * VARIANCE_THRESHOLD)
        M = [MixtureModel(kernels) for i in 1:n, j in 1:m]
        new(M, MixtureParams(), cut)
    end
end

function apply!(M::MarkovModel, image)
    update!(M.cut, bg_energy!.(M.background, M.params, image))
    increment_frames!(M.params)
    segment!(M.cut)
end

function apply_mrf!(cut, background, params, image, smoothing=2.0)
    energy = bg_energy!.(background, params, image)
    increment_frames!(params)
    reset!(cut, energy,
        VARIANCE_THRESHOLD,
        smoothing * VARIANCE_THRESHOLD
    )
    segment!(cut)
end
