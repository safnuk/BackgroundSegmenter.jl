using BackgroundSegmenter
using Base.Test

K = Gaussian(0.0, 1.0)
@test neglog_likelihood(0.0, K) ≈ 0 
@test neglog_likelihood(0, K) ≈ 0 

L = Gaussian(0, 2)
@test neglog_likelihood(1, L) ≈ 0.5

M = Gaussian(0, 4, 4)
@test M.sortkey == 2

model = MixtureModel([K, L, M])
@test model[1].σ² == M.σ²
model[end] = L
@test model[3].σ² == L.σ²

model = MixtureModel(5)
@test size(model) == (5,)
@test model[1].μ == 0 

model = MixtureModel([Gaussian(0, 4, 1), Gaussian(30, 4, 0.5), Gaussian(70, 4, 0.48), Gaussian(120, 4, .02)])
out = apply!(model, 0)
@test out == 0
@test model[1].ω > 0.5
@test model[2].ω < 0.25
@test model[3].ω < 0.24

model = MixtureModel([Gaussian(0, 4, 1), Gaussian(30, 4, 0.5), Gaussian(70, 4, 0.48), Gaussian(120, 4, .02)])
model.number_frames = 200
out = apply!(model, 120)
@test out == 1
@test model[1].ω < 0.5
@test model[2].ω < 0.25
@test model[3].ω < 0.24
@test model[4].ω > 0.01
@test bg_energy!(model, 30) == 0
@test bg_energy!(model, 120) ≈ neglog_likelihood(120, model[3])

bg_input = [0x90 0x30 ; 0x70 0x120]
fg_input = [0x30 0x70; 0x120 0x20]
output = [0x0 0x0; 0x0 0x0]
model = [MixtureModel(3) for i in 1:2, j in 1:2]
apply!.(model, bg_input)
for i in 1:BackgroundSegmenter.INITIALIZATION_WINDOW
    output .+= apply!.(model, bg_input)
end
@test output == [0x0 0x0; 0x0 0x0]
output .= apply!.(model, fg_input)
@test output == [0x1 0x1; 0x1 0x1]
