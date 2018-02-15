__precompile__()

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
    @time fgbg[i, :, :] = label_components(apply_mrf!(cut, mm, (@view V[i, :, :]), 2.0), 8)
end

V = load("data/lamprey2.npy");
V = V[:, 300:450, 150:400];
fgbg = zeros(V);
(t, n, m) = size(V)
M = [MixtureModel(5) for i in 1:n, j in 1:m];
fgbg = zeros(V);
fgbg[1, :, :] = label_components(apply!.(M, @view V[1, :, :]), 10);
@time fgbg[2, :, :] = apply!.(M, @view V[2, :, :]);
@time fgbg[2, :, :] = label_components((@view fgbg[2, :, :]), 10);
@time for i in 3:t-1
    fgbg[i, :, :] = label_components(apply!.(M, @view V[i, :, :]), 10);
end
@time fgbg[t, :, :] = apply!.(M, @view V[t, :, :]);
@time fgbg[t, :, :] = label_components((@view fgbg[t, :, :]), 10);
B = permutedims(V .* fgbg, [2, 3, 1]);

mm =  [MixtureModel(5) for i in 1:n, j in 1:m];
cut = MinCut(n, m);
 for i in 1:5
    time_results(V, fgbg, cut,  mm, i)
end
@time  for i in 6:t-1
    fgbg[i, :, :] = label_components(apply_mrf!(cut, mm, (@view V[i, :, :]), 2.0), 12)
end
time_results(V, fgbg, cut,  mm, t)

A = permutedims(V .* fgbg, [2, 3, 1])
C = zeros(2n, 2m, t);
C[1:n, 1:m, :] = A;
C[n+1:end, 1:m, :] = B;
C[1:n, m+1:end, :] = permutedims(V, [2, 3, 1])

#using ImageView, GtkReactive, Colors
#imshow(C)
