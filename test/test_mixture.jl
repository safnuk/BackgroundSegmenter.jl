using BackgroundSegmenter
using Base.Test

K = Gaussian(0.0, 1.0)
@test neglog_likelihood(K, 0.0) ≈ 0
@test neglog_likelihood(K, 0) ≈ 0

L = Gaussian(0, 2)
@test neglog_likelihood(L, 1) ≈ 0.5

M = Gaussian(0, 4, 4)
@test M.sortkey == 2

bg_input = [0x90 0x10 ; 0x50 0x120]
fg_input = [0x30 0x70; 0x120 0x20]
output = [0x0 0x0; 0x0 0x0]
model = MixtureModel(3, 2, 2)
apply!(model, bg_input)
for i in 1:BackgroundSegmenter.INITIALIZATION_WINDOW
    output .+= apply!(model, bg_input)
end
@test output == [0x0 0x0; 0x0 0x0]
output .= apply!(model, fg_input)
@test output == [0x0 0x1; 0x1 0x0]
