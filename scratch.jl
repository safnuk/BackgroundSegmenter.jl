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
        
function fgbg_markov(infile, outfile, k=5)
    V = np.load(infile)
    (t, n, m) = size(V)
    mm = MarkovModel(n, m, k, 2.0)
    fgbg = zeros(V)
    for i in 1:t
        fgbg[i, :, :] = apply!(mm, V[i, :, :])
    end
    np.save(outfile, fgbg)
end
