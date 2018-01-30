using Base.Test

using BackgroundSegmenter

A = [1 0 1 1 ; 1 0 1 0 ; 0 1 1 0; 1 1 0 0]
out = label_components(A, 3)
@test out == [0 0 1 1 ; 0 0 1 0 ; 0 1 1 0; 1 1 0 0]
out = label_components(A, 2)
@test out == A
out = label_components(A, 8)
@test out == zeros(A)
