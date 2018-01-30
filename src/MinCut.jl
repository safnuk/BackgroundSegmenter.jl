using DataStructures
using IterTools

GridGraph = Array{Float64, 3}
Node = Tuple{Int64, Int64}
Edge = Tuple{Node, Node}
Path = Vector{Node}

const _SOURCE = -2 
const SOURCE = (_SOURCE, _SOURCE)
const _SINK = -1
const SINK = (_SINK, _SINK)
const _NULL = 0
const NULL = (_NULL, _NULL)

struct MinCut
    orphans::Set
    active::Queue
    tree::DefaultDict{Node, Int}
    parent::DefaultDict{Node, Node}

    function MinCut() 
        orphans = Set{Node}()
        active = Queue(Node)
        enqueue!(active, SOURCE)
        enqueue!(active, SINK)
        tree = DefaultDict{Node, Int}(_NULL)
        parent = DefaultDict{Node, Node}(NULL)
        tree[SOURCE] = _SOURCE
        tree[SINK] = _SINK
        new(orphans, active, tree, parent)
    end
end

function segment(G::GridGraph)
    cut = mincut(G)
    (n, m, k) = size(G)
    fgbg = zeros(UInt8, n, m)
    for i in 1:n, j in 1:m
        if cut.tree[(i,j)] == _SOURCE
            fgbg[i, j] = zero(UInt8)
        else
            fgbg[i, j] = one(UInt8)
        end
    end
    return fgbg
end

function mincut(G::GridGraph)
    cut = MinCut()
    while true
    #for n in 1:5
        path = grow!(cut, G)
        if length(path) == 0
            return cut
        end
        augment!(cut, G, path)
        adopt!(cut, G)
    end
end

function grow!(cut::MinCut, G::GridGraph)
    while length(cut.active) > 0
        p = front(cut.active)
        for q in neighbors(G, p)
            if tree_capacity(cut, G, p, q) < eps()
                continue
            elseif cut.tree[q] == _NULL
                cut.tree[q] = cut.tree[p]
                cut.parent[q] = p
                enqueue!(cut.active, q)
            elseif cut.tree[q] != cut.tree[p]
                # should we dequeue p before this?
                return trace_path(cut, p, q)
            end
        end
        dequeue!(cut.active)
    end
    return Path(0)
end

function neighbors(G::GridGraph, p::Node)
    (a, b) = p
    (n, m) = size(G)
    if p == SINK || p == SOURCE
        return [(i, j) for i in 1:n, j in 1:m]
    else
        possible = [(a+1, b), (a-1, b), (a, b-1), (a, b+1)]
        filter!(x -> 1<=x[1] && x[1]<=n && 1<=x[2] && x[2]<=m, possible)
        return [possible; SINK; SOURCE]
    end
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

function tree_capacity(cut::MinCut, G::GridGraph, p::Node, q::Node)
    if cut.tree[p] == _SOURCE
        return edge_weight(G::GridGraph, p, q)
    else
        return edge_weight(G::GridGraph, q, p)
    end
end

function edge_weight(G::GridGraph, p::Node, q::Node)
    (n, m, k) = weight_index(p, q)
    if k < 0
        return 0
    else
        return G[n, m, k]
    end
end

function weight_index(p::Node, q::Node)
    (x, y) = p
    (z, w) = q
    if p == SINK
        return (z, w, -1)
    elseif q == SOURCE
        return (x, y, -1)
    elseif p == SOURCE
        return (z, w, 6)
    elseif q == SINK
        idx = 5 
    elseif x < z
        idx = 4
    elseif x > z
        idx = 2
    elseif y < w
        idx = 3
    elseif y > w
        idx = 1
    end
    return (x, y, idx)
end

function augment!(cut::MinCut, G::GridGraph, path::Path)
    Δ = bottleneck_capacity(G, path)
    saturated_edges = update_residual_capacity!(G, path, Δ)
    for edge in saturated_edges
        (p, q) = edge
        if cut.tree[p] != cut.tree[q] || cut.tree[p] == _NULL
            continue
        elseif cut.tree[p] == _SOURCE == cut.tree[q]
            cut.parent[q] = NULL
            push!(cut.orphans, q)
        elseif cut.tree[p] == _SINK == cut.tree[q]
            cut.parent[p] = NULL
            push!(cut.orphans, p)
        end
    end
end

function bottleneck_capacity(G::GridGraph, path::Path)
    state = start(path)
    (p, state) = next(path, state)
    capacity = Inf
    while !done(path, state)
        (q, state) = next(path, state)
        capacity = min(capacity, edge_weight(G, p, q))
        p = q
    end
    return capacity
end

function update_residual_capacity!(G::GridGraph, path::Path, Δ::Float64)
    saturated = Queue(Edge)
    state = start(path)
    (p, state) = next(path, state)
    while !done(path, state)
        (q, state) = next(path, state)
        if update_residual_capacity!(G, p, q, Δ) < eps()
            enqueue!(saturated, (p, q))
        end
        p = q
    end
    return saturated
end

function update_residual_capacity!(G::GridGraph, p::Node, q::Node, Δ::Float64)
    (n, m, k) = weight_index(q, p)
    if k > 0
        G[n, m, k] += Δ
    end
    (n, m , k) = weight_index(p, q)
    if k > 0
        G[n, m ,k] -= Δ
        return G[n, m, k]
    else
        return 0.0
    end
end

function adopt!(cut::MinCut, G::GridGraph)
    while !isempty(cut.orphans)
        p = pop!(cut.orphans)
        attempt_tree_graft!(cut, G, p)
    end
end

function attempt_tree_graft!(cut::MinCut, G::GridGraph, p::Node)
    for q in neighbors(G, p)
        if isvalid_parent_of(cut, G, q, p)
            cut.parent[p] = q
            return
        end
    end
    for q in filter(x-> cut.tree[x] == cut.tree[p], neighbors(G, p))
        if tree_capacity(cut, G, q, p) > 0
            enqueue!(cut.active, q)
        end
        if cut.parent[q] == p
            push!(cut.orphans, q)
            cut.parent[q] = NULL
        end
    end
    cut.tree[p] = _NULL
    remove!(cut.active, p)
end

function remove!(q::Queue{T}, node::T) where T
    p = Queue(T)
    while !isempty(q)
        r = dequeue!(q)
        if r != node
            enqueue!(p, r)
        end
    end
    for r in p
        enqueue!(q, r)
    end
    return q
end

function isvalid_parent_of(cut::MinCut, G::GridGraph, q::Node, p::Node)
    if cut.tree[q] != cut.tree[p]
        return false
    elseif tree_capacity(cut, G, q, p) < eps()
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
