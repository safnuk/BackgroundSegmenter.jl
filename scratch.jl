__precompile__()

using PyCall
using ImageView, GtkReactive, Colors, Images

using BackgroundSegmenter

@pyimport skvideo as skv
@pyimport numpy as np
@pyimport skvideo.io as skvio

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

function time_results(V, fgbg, M, t)
    @time for i in 2:t
        apply!(M, (@view V[:, :, i]), (@view fgbg[:, :, i]));
    end
end

# V = permutedims(load("data/lamprey2.npy"), [2, 3, 1]);
# V = V[300:450, 150:400, :];
V = permutedims(skvio.vread("data/lamprey4.avi")[:, :, :, 3], [2, 3, 1]);
f1 = BackgroundSegmenter.Filter(1, 0.2);
f2 = BackgroundSegmenter.Filter(3, 0.08);
f3 = BackgroundSegmenter.Filter(5, 0.04);
fgbg = zeros(V);
(n, m, t) = size(V)
M = MixtureModel(n, m);
apply!(M, (@view V[:, :, 1]), (@view fgbg[:, :, 1]));
time_results(V, fgbg, M, t)

for i in 1:t
    fgbg[i, :, :] .= filter_components(fgbg[i, :, :], 12)
end
@time fgbg .= filter_components(fgbg, 40);
@time fgbg1 = apply(f1, fgbg);
A0 = V;
A1 = scale(V, fgbg1);
guidict = imshow(hcat(A0, A1))
# @time fgbg2 = apply(f2, fgbg);
# @time fgbg3 = apply(f3, fgbg);
# A2 = scale(V, fgbg2);
# A3 = scale(V, fgbg3);

# C = zeros(UInt8, 2n, 2m, t);
# C[1:n, 1:m, :] = A0;
# C[n+1:end, 1:m, :] = A2;
# C[1:n, m+1:end, :] = A1;
# C[n+1:end, m+1:end, :] = A3;

# idx = annotate!(guidict, AnnotationBox((20, 30), (40, 20), linewidth=2, color=RGB(0, 1, 0)))
# mm =  [MixtureModel(5) for i in 1:n, j in 1:m];
# cut = MinCut(n, m);
# for i in 1:5
#     time_results(V, fgbg, cut,  mm, i)
# end
# @time  for i in 6:t-1
#     fgbg[i, :, :] = label_components(apply_mrf!(cut, mm, (@view V[i, :, :]), 2.0), 12)
# end
# time_results(V, fgbg, cut,  mm, t)

# A = permutedims(V .* fgbg, [2, 3, 1])
