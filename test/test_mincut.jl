grid = ones(Float64, 4, 4)
grid[1:2, :] = 4.0
grid[3:4, :] = 0.2
cut = MinCut(size(grid)...)
reset!(cut, grid, 1.0, 1.0)
out = ones(UInt8, 4, 4)
out[3:4, :] = zero(UInt8)
@test segment!(cut) == out
