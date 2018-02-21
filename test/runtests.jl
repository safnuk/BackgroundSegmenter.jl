using Base.Test
using BackgroundSegmenter

tic()
@time @testset "Test queue" begin include("test_queue.jl") end
@time @testset "Test mixture model" begin include("test_mixture.jl") end
@time @testset "Test component finder" begin include("test_components.jl") end
@time @testset "Test MinCut algorithm" begin include("test_mincut.jl") end
@time @testset "Test Path data structure" begin include("test_path.jl") end
toc()
