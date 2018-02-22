using Base.Test

using BackgroundSegmenter

A = [1 0 1 1; 1 0 1 0; 0 1 1 0; 1 1 0 0]
B = [1 0 1 1; 0 0 1 0; 0 0 0 0; 0 0 1 1]
C = zeros(Int64, 2, 4, 4)
C[1, :, :] = A
C[2, :, :] = B
out = filter_components(A, 3)
@test out == [0 0 1 1 ; 0 0 1 0 ; 0 1 1 0; 1 1 0 0]
out = filter_components(A, 2)
@test out == [1 0 2 2 ; 1 0 2 0 ; 0 2 2 0; 2 2 0 0]
out = filter_components(A, 8)
@test out == zeros(A)
out = filter_components(C, 2)
@test sum(out) == 29
out = filter_components(C, 3)
@test sum(out) == 23
out = filter_components(C, 4)
@test sum(out) == sum(C) - 5
out = filter_components(C, 12)
@test out == zeros(C)
