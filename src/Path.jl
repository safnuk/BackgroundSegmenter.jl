using DataStructures
import DataStructures: enqueue!, reverse_iter, back
import Base: length, next, start, done

struct HalfPath
    q::Queue{Int}
    HalfPath() = new(Queue(Int))
end

@DataStructures.delegate HalfPath.q [length, back, enqueue!, next, start, done, reverse_iter]

struct Path
    head::HalfPath
    tail::HalfPath
end
empty_path() = Path(HalfPath(), HalfPath())

Base.length(P::Path) = length(P.head) + length(P.tail)
Base.start(P::Path) = (false, start(reverse_iter(P.head)))
Base.done(P::Path, state) = state[1] && done(P.tail, state[2])
function Base.next(P::Path, state)
    in_tail, s = state
    if in_tail
        item, ns = next(P.tail, s)
    elseif done(reverse_iter(P.head), s)
        ns = start(P.tail)
        item, ns = next(P.tail, ns)
        in_tail = true
    else
        item, ns = next(reverse_iter(P.head), s)
    end
    return (item, (in_tail, ns))
end
Base.eltype(::Type{Path}) = Int
