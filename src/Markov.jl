struct MarkovModel
    background :: Array{MixtureModel, 2}
    clique_potential :: Float64
    function MarkovModel(n::Int, m::Int, kernels=5, clique_potential=2.0)
        M = [MixtureModel(kernels) for i in 1:n, j in 1:m]
        new(M, clique_potential)
    end
end

function apply!(M::MarkovModel, image)
    (n, m) = size(image)
    G = ones(Float64, n, m, 6)
    G[:, :, 1:4] = M.clique_potential * INITIAL_VARIANCE 
    G[:, :, 6] = INITIAL_VARIANCE
    G[:, :, 5] = bg_energy!.(M.background, image)
    segment(G)
end
