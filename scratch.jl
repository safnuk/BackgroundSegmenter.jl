using BackgroundSegmenter
using PyCall

@pyimport numpy as np

function fgbg_mog(infile, outfile, k=5)
    V = np.load(infile)
    (t, n, m) = size(V)
    M = [MixtureModel(k) for i in 1:n, j in 1:m]
    fgbg = zeros(V)
    for i in 1:t
        fgbg[i, :, :] = label_components(apply!.(M, V[i, :, :]), 10)
    end
    np.save(outfile, fgbg)
end
        
function load(infile, limit=0)
    V = np.load(infile)
    (t, n, m) = size(V)
    if limit == 0
        return V
    end
    left = convert(Int, round((n - limit) / 2))
    top = convert(Int, round((m - limit) / 2))
    return V[:, left:(left+limit), top:(top+limit)]
end

function fgbg_markov(infile, outfile, k=5, limit=0)
    V = np.load(infile)
    (t, n, m) = size(V)

    mm = MarkovModel(n, m, k, 2.0)
    fgbg = zeros(V)
    for i in 1:t
        fgbg[i, :, :] = apply!(mm, @view V[i, :, :])
    end
    np.save(outfile, fgbg)
end

V = load("data/lamprey1.npy", 20);
mm = MarkovModel(V[1,:,:], 5);
fgbg = zeros(V);
@time fgbg[1, :, :] = apply!(mm, @view V[1, :, :])
for i in 1:t
    fgbg[i, :, :] = apply!(mm, @view V[i, :, :])
end
