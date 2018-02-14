__precompile__()

using ImageView
using PyCall

using BackgroundSegmenter

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

    # mm = MarkovModel(n, m, k, 2.0)
    M = [MixtureModel(k) for i in 1:n, j in 1:m]
    fgbg = zeros(V)
    for i in 1:t
        fgbg[i, :, :] = apply!(mm, @view V[i, :, :])
    end
    np.save(outfile, fgbg)
end


function time_results(V, fgbg, cut, mm, i)
    @time fgbg[i, :, :] = label_components(apply_mrf!(cut, mm, V[i, :, :], 4.0), 8)
end

V = load("data/lamprey2.npy");
V = V[1:10, 300:450, 150:400]
fgbg = zeros(V);
(t, n, m) = size(V)
mm =  [MixtureModel(5) for i in 1:n, j in 1:m];
cut = MinCut(n, m);

for i in 1:10
    time_results(V, fgbg, cut,  mm, i)
end

A = V .* fgbg
imshow(A)
