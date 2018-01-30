G = ones(Float64, 4, 4, 6)
G[1:2, :, 5] = 0.2
G[3:4, :, 5] = 4.0
out = ones(UInt8, 4, 4)
out[1:2, :] = zero(UInt8)
segment(G)
@test segment(G) == out
