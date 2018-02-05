using DataStructures
using IterTools
import Base.length
import Base.size

GridGraph = Array{Float64, 2}

Node = Int64
Edge = Tuple{Node, Node}
Path = Vector{Node}

const SOURCE = 1 
const SINK = 2
const NULL = 0

struct MinCut
    orphans::Set{Node}
    active::OrderedSet{Node}
    tree::Vector{Node}
    parent::Vector{Node}
    edges::Array{Node, 2}
    weights::Array{Float64, 2}
    _n::Int
    _m::Int

    function MinCut(G::GridGraph, sink, node) 
        num_nodes = length(G) + 2
        orphans = Set{Node}()
        active = OrderedSet{Node}([SINK])
        for k in 3:num_nodes
            push!(active, k)
        end
        tree = Vector{Node}(num_nodes)
        tree[:] = SOURCE
        tree[SINK] = SINK
        parent = Vector{Node}(num_nodes)
        parent[:] = SOURCE
        parent[SOURCE] = NULL
        parent[SINK] = NULL
        (edges, weights) = build_edges(G, sink, node)
        (n, m) = size(G)
        new(orphans, active, tree, parent, edges, weights, n, m)
    end
end

size(cut::MinCut) = (cut._n, cut._m)
length(cut::MinCut) = length(cut.tree)

function build_edges(G::GridGraph, sink, node)
    offsets = [(-1, 0), (0, -1), (1, 0), (0,1)]
    (n, m) = size(G)
    num_nodes = n * m + 2
    weights = Array{Float64, 2}(num_nodes, 6)
    edges = zeros(Node, num_nodes, 6)
    weights[:, 3:6] = node
    weights[3:end, SOURCE] = reshape(G, n*m)
    weights[3:end, SINK] = sink
    weights[1:2, 1:2] = 0.0
    edges[:, SOURCE] = SOURCE
    edges[:, SINK] = SINK
    for i in 1:n, j in 1:m
        k = linearize((i, j), n)
        multi_idx = [(i, j) .+ x for x in offsets]
        multi_idx = [(r, s) for (r, s) in multi_idx if
                     r >= 1 && s >= 1 && r <= n && s <= m]
        idx = linearize.(multi_idx, n)
        edges[k, 3:length(idx)+2] = idx
        weights[k, length(idx)+3:end] = 0.0
        for (q, (r, s)) in enumerate(multi_idx)
            if r < 1 || s < 1 || r > n || s > m
                weights[k, q+2] = 0
            end
        end
    end
    return (edges, weights)
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
        if cut.tree[i+2] == SOURCE
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
    while length(cut.active) > 0
        # get the first item in the queue
        state = start(cut.active)
        (p, state) = next(cut.active, state)
        for q in viable_neighbors(cut, p)
            if cut.tree[q] == NULL
                cut.tree[q] = cut.tree[p]
                cut.parent[q] = p
                push!(cut.active, q)
            elseif cut.tree[q] != cut.tree[p]
                return trace_path(cut, p, q)
            end
        end
        delete!(cut.active, p)
    end
    return Path(0)
end

function viable_neighbors(cut::MinCut, p::Node)
    if p == SOURCE
        # TODO: memoize this process
        return filter(x -> cut.weights[x, SOURCE] > eps(), 3:length(cut))
    elseif p == SINK
        # TODO: memoize this process
        return filter(x -> cut.weights[x, SINK] > eps(), 3:length(cut))
    elseif cut.tree[p] == SOURCE
        return [cut.edges[p, x] for x in [SINK; 3:6] if cut.weights[p, x] > eps()]
    else
        return [cut.edges[p, x] for x in [SOURCE; 3:6] if
                cut.edges[p, x] > 0 && capacity(cut, cut.edges[p, x], p) > eps()]
    end
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
    q_path = trace_to_root(cut, q)
    p_path = trace_to_root(cut, p)
    if back(q_path) == SOURCE
        return [z for z in chain(reverse_iter(q_path), p_path)]
    else
        return [z for z in chain(reverse_iter(p_path), q_path)]
    end
end

function trace_to_root(cut::MinCut, p::Node)
    path = Queue(Node)
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
    saturated_edges = update_residual_capacity!(cut, path, Δ)
    for edge in saturated_edges
        (p, q) = edge
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
    @assert cap > 0
    return cap
end

function update_residual_capacity!(cut::MinCut, path::Path, Δ::Float64)
    saturated = Queue(Edge)
    state = start(path)
    (p, state) = next(path, state)
    while !done(path, state)
        (q, state) = next(path, state)
        if update_residual_capacity!(cut, p, q, Δ) < eps()
            enqueue!(saturated, (p, q))
        end
        p = q
    end
    return saturated
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
            push!(cut.active, q)
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
        orig = origin(cut, q)
        return orig == SOURCE || orig == SINK
    end
end

function origin(cut::MinCut, p::Node)
    while cut.parent[p] != NULL
        p = cut.parent[p]
    end
    return p
end
