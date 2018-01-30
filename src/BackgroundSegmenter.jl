module BackgroundSegmenter
export
    apply!,
    bg_energy!,
    Gaussian,
    label_components,
    MarkovModel,
    MixtureModel,
    mincut,
    neglog_likelihood,
    segment

include("Components.jl")
include("MinCut.jl")
include("Mixture.jl")
include("Markov.jl")
end # module
