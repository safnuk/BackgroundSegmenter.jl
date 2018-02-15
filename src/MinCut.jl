using DataStructures
using IterTools
import Base.length
import Base.size
import DataStructures.reset!

GridGraph = Array{Float64, 2}

Node = Int64
Edge = Tuple{Node, Node}

const SOURCE = 1
const SINK = 2
const NULL = 0

mutable struct MinCut
    orphans::Set{Node}
    active::ActiveQueue
    tree::Vector{Node}
    parent::Vector{Node}
    edges::Array{Node, 2}
    weights::Array{Float64, 2}
    timestamps::Vector{Int}
    distances::Vector{Int}
    current_time::Int
    _n::Int
    _m::Int
    _p1::HalfPath
    _p2::HalfPath
    saturated::HalfPath

    function MinCut(n, m)
        num_nodes = n * m + 2
        orphans = Set{Node}()
        active = ActiveQueue(num_nodes)
        tree = Vector{Node}(num_nodes)
        parent = Vector{Node}(num_nodes)
        weights = Array{Float64, 2}(num_nodes, 6)
        edges = build_edges(n, m)
        times = zeros(Int, num_nodes)
        distances = zeros(Int, num_nodes)
        new(orphans, active, tree, parent, edges, weights, times, distances, 0, n, m, HalfPath(), HalfPath(), HalfPath())
    end
end

size(cut::MinCut) = (cut._n, cut._m)
length(cut::MinCut) = length(cut.tree)

function reset!(cut::MinCut, G::GridGraph, sink, node)
    @assert (cut._n, cut._m) == size(G)
    clear!(cut.active)
    enqueue!(cut.active, SINK)
    cut.tree[SOURCE] = SOURCE
    cut.tree[SINK] = SINK
    cut.parent[SOURCE] = NULL
    cut.parent[SINK] = NULL
    cut.weights[:, 3:6] = node
    cut.weights[3:end, SOURCE] = reshape(G, cut._n*cut._m)
    cut.weights[3:end, SINK] = sink
    cut.weights[1:2, 1:2] = 0.0
    cut.timestamps[:] = 1
    for k in 3:length(cut)
        if cut.weights[k, SOURCE] < cut.weights[k, SINK] - eps()
            cut.weights[k, SINK] -= cut.weights[k, SOURCE]
            cut.weights[k, SOURCE] = 0.0
            cut.tree[k] = SINK
            cut.parent[k] = SINK
            cut.distances[k] = 1
        elseif cut.weights[k, SOURCE] > cut.weights[k, SINK] + eps()
            cut.weights[k, SOURCE] -= cut.weights[k, SINK]
            cut.weights[k, SINK] = 0.0
            cut.tree[k] = SOURCE
            cut.parent[k] = SOURCE
            enqueue!(cut.active, k)
            cut.distances[k] = 1
        else
            cut.weights[k, SOURCE] = 0.0
            cut.weights[k, SINK] = 0.0
            cut.tree[k] = NULL
            cut.parent[k] = NULL
            cut.distances[k] = -1
        end
    end
    return cut
end

function build_edges(n, m)
    offsets = [(-1, 0), (0, -1), (1, 0), (0,1)]
    num_nodes = n * m + 2
    edges = zeros(Node, num_nodes, 6)
    edges[:, SOURCE] = SOURCE
    edges[:, SINK] = SINK
    for i in 1:n, j in 1:m
        k = linearize((i, j), n)
        multi_idx = [(i, j) .+ x for x in offsets]
        multi_idx = [(r, s) for (r, s) in multi_idx if
                     r >= 1 && s >= 1 && r <= n && s <= m]
        idx = linearize.(multi_idx, n)
        edges[k, 3:length(idx)+2] = idx
        #weights[k, length(idx)+3:end] = 0.0
        #for (q, (r, s)) in enumerate(multi_idx)
        #    if r < 1 || s < 1 || r > n || s > m
        #        weights[k, q+2] = 0
        #    end
        #end
    end
    return edges
end

linearize(idx, n) = idx[1] + (idx[2]-1)*n + 2

function update!(cut::MinCut, source::GridGraph)
    cut.weights[3:end, SOURCE] = reshape(source, length(source))
end

function segment!(cut::MinCut)
    mincut!(cut)
    (n, m) = size(cut)
    fgbg = zeros(UInt8, n, m)
    for i in 1:length(fgbg)
        if cut.tree[i+2] == SINK
            fgbg[i] = zero(UInt8)
        else
            fgbg[i] = one(UInt8)
        end
    end
    return fgbg
end

function mincut!(cut::MinCut)
    while true
        path = grow!(cut)
        if length(path) == 0
            return cut
        end
        augment!(cut, path)
        adopt!(cut)
    end
end

function grow!(cut::MinCut)
    while !isempty(cut.active)
        p = front(cut.active)
        for q in neighbors(cut, p)
            if cut.tree[p] == SOURCE && q == SOURCE
                continue
            elseif cut.tree[p] == SINK && q == SINK
                continue
            elseif tree_capacity(cut, p, q) < eps()
                continue
            elseif cut.tree[q] == NULL
                cut.tree[q] = cut.tree[p]
                cut.parent[q] = p
                enqueue!(cut.active, q)
            elseif cut.tree[q] != cut.tree[p]
                return trace_path(cut, p, q)
            end
        end
        dequeue!(cut.active)
    end
    return empty_path()
end

function neighbors(cut::MinCut, p::Node)
    if p == SOURCE || p == SINK
        return 3:length(cut)
    else
        m = 0
        for m in 3:6
            if cut.edges[p, m] == 0
                return @view cut.edges[p, 1:m-1]
            end
        end
        return @view cut.edges[p, :]
    end
end

function tree_capacity(cut::MinCut, p::Node, q::Node)
    if cut.tree[p] == SOURCE
        return capacity(cut, p, q)
    else
        return capacity(cut, q, p)
    end

end

function capacity(cut::MinCut, p::Node, q::Node)
    if p == SOURCE
        idx = SOURCE
        p = q
    elseif q == SINK
        idx = SINK
    else
        idx = edge_index(cut, p, q)
    end
    return cut.weights[p, idx]
end

function edge_index(cut::MinCut, p::Node, q::Node)
    for (n, x) in enumerate(@view cut.edges[p, 3:end])
        if x == q
            return n + 2
        end
    end
    @assert false
end

function trace_path(cut::MinCut, p::Node, q::Node)
    empty!(cut._p1)
    empty!(cut._p2)
    q_path = trace_to_root!(cut._p1, cut, q)
    p_path = trace_to_root!(cut._p2, cut, p)
    if back(q_path) == SOURCE
        return Path(q_path, p_path)
    else
        return Path(p_path, q_path)
    end
end

function trace_to_root!(path::HalfPath, cut::MinCut, p::Node)
    enqueue!(path, p)
    parent = cut.parent[p]
    while parent != NULL
        enqueue!(path, parent)
        parent = cut.parent[parent]
    end
    return path
end

function augment!(cut::MinCut, path::Path)
    Δ = bottleneck_capacity(cut, path)
    update_residual_capacity!(cut, path, Δ)
    while !isempty(cut.saturated)
        p = dequeue!(cut.saturated)
        q = dequeue!(cut.saturated)
        if cut.tree[p] != cut.tree[q] || cut.tree[p] == NULL
            continue
        elseif cut.tree[p] == SOURCE == cut.tree[q]
            cut.parent[q] = NULL
            push!(cut.orphans, q)
        elseif cut.tree[p] == SINK == cut.tree[q]
            cut.parent[p] = NULL
            push!(cut.orphans, p)
        end
    end
end

function bottleneck_capacity(cut::MinCut, path::Path)
    state = start(path)
    (p, state) = next(path, state)
    cap = Inf
    while !done(path, state)
        (q, state) = next(path, state)
        cap = min(cap, capacity(cut, p, q))
        p = q
    end
    return cap
end

function update_residual_capacity!(cut::MinCut, path::Path, Δ::Float64)
    state = start(path)
    (p, state) = next(path, state)
    while !done(path, state)
        (q, state) = next(path, state)
        if update_residual_capacity!(cut, p, q, Δ) < eps()
            push!(cut.saturated, p, q)
        end
        p = q
    end
end

function update_residual_capacity!(cut::MinCut, p::Node, q::Node, Δ::Float64)
    if p == SOURCE
        cut.weights[q, SOURCE] -= Δ
        return cut.weights[q, SOURCE]
    elseif q == SINK
        cut.weights[p, SINK] -= Δ
        return cut.weights[p, SINK]
    end
    idx = edge_index(cut, p, q)
    opp_idx = edge_index(cut, q, p)
    cut.weights[q, opp_idx] += Δ
    cut.weights[p, idx] -= Δ
end

function adopt!(cut::MinCut)
    cut.current_time += 1
    while !isempty(cut.orphans)
        p = pop!(cut.orphans)
        attempt_tree_graft!(cut, p)
    end
end

function attempt_tree_graft!(cut::MinCut, p::Node)
    neighborhood = neighbors(cut, p)
    for q in neighborhood
        if isvalid_parent_of(cut, q, p)
            cut.parent[p] = q
            return
        end
    end
    for q in filter(x-> cut.tree[x] == cut.tree[p], neighborhood)
        if has_capacity(cut, q, p)
            enqueue!(cut.active, q)
        end
        if cut.parent[q] == p
            push!(cut.orphans, q)
            cut.parent[q] = NULL
        end
    end
    cut.tree[p] = NULL
    delete!(cut.active, p)

end

function has_capacity(cut::MinCut, p::Node, q::Node)
    if cut.tree[p] == SOURCE
        return capacity(cut, p, q) > eps()
    else
        return capacity(cut, q, p) > eps()
    end
end

function isvalid_parent_of(cut::MinCut, q::Node, p::Node)
    if cut.tree[q] != cut.tree[p]
        return false
    elseif !has_capacity(cut, q, p)
        return false
    else
        d = distance_to_origin!(cut, q)
        return d >= 0
    end
end

function distance_to_origin!(cut::MinCut, p::Node)
    # TODO: Eliminate recursion?
    if cut.timestamps[p] == cut.current_time
        return cut.distances[p]
    end
    cut.timestamps[p] = cut.current_time
    if p == SOURCE || p == SINK
        return 0
    elseif cut.parent[p] == NULL
        cut.distances[p] = -1
        return -1
    else
        d = distance_to_origin!(cut, cut.parent[p])
        if d < 0
            cut.distances[p] = d
        else
            cut.distances[p] = d + 1
        end
    end
    return cut.distances[p]
end
