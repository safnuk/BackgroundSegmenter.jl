using DataStructures

import Base: isempty, length, delete!
import DataStructures: enqueue!, dequeue!, front


mutable struct ActiveQueue
    _v::Vector{Int}
    removed::Vector{Int}
    active::Vector{Bool}
    head::Int
    tail::Int
    prior_state::Tuple{Int, Int}
    ActiveQueue(max_index) = new(
        Vector{Int}(2max_index), zeros(max_index), falses(max_index), 1, 0, (-1, 0))
end

function clear!(q::ActiveQueue)
    for n in 1:length(q.removed)
        q.removed[n] = 0
        q.active[n] = false
    end
    q.head = 1
    q.tail = 0
    q.prior_state = (-1, 0)
end

function enqueue!(q::ActiveQueue, x::Int)
    # Warning: Does not check bounds of x
    if q.active[x] # don't add an existing item to queue
        return q
    elseif q.tail < length(q._v)
        q.active[x] = true
        q.tail += 1
        q._v[q.tail] = x
    else
        rebase!(q)
        enqueue!(q, x)
    end
end

function dequeue!(q::ActiveQueue)
    result = front(q)
    q.head += 1
    q.active[result] = false
    return result
end

function rebase!(q::ActiveQueue)
    println("Rebasing")
    current_index = 1
    for idx in q.head:q.tail
        x = q._v[idx]
        if q.removed[x] > 0
            q.removed[x] -= 1
        else
            q._v[current_index] = x
            current_index += 1
        end
    end
    q.head = 1
    q.tail = current_index
end

function isempty(q::ActiveQueue)
    while q.head <= q.tail && q.removed[q._v[q.head]] > 0
        q.removed[q._v[q.head]] -= 1
        q.head += 1
    end
    return q.head > q.tail
end

function front(q::ActiveQueue)
    if isempty(q)
        throw(ArgumentError("ActiveQueue must be non-empty"))
    end
    return q._v[q.head]
end

function delete!(q::ActiveQueue, x::Int)
    if q.active[x]
        q.active[x] = false
        q.removed[x] += 1
    end
end
