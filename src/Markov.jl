struct MarkovModel
    background :: Array{MixtureModel, 2}
    cut::MinCut
    function MarkovModel(G, kernels=5, smoothing=2.0)
        (n, m) = size(G)
        cut = MinCut(zeros(Float64, n, m),
                     VARIANCE_THRESHOLD,
                     smoothing * VARIANCE_THRESHOLD)
        M = [MixtureModel(kernels) for i in 1:n, j in 1:m]
        new(M, cut)
    end
end

function apply!(M::MarkovModel, image)
    update!(M.cut, bg_energy!.(M.background, image))
    segment!(M.cut)
end

function apply_mrf!(cut, background, image, smoothing=2.0)
    energy = bg_energy!.(background, image)
    reset!(cut, energy,
        VARIANCE_THRESHOLD,
        smoothing * VARIANCE_THRESHOLD
    )
    segment!(cut)
end
