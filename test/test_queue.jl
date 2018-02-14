q = ActiveQueue(6)
enqueue!(q, 1)
@test front(q) == 1
enqueue!(q, 2)
enqueue!(q, 3)
enqueue!(q, 4)
@test dequeue!(q) == 1
@test dequeue!(q) == 2
delete!(q, 4)
delete!(q, 3)
enqueue!(q, 5)
enqueue!(q, 6)
enqueue!(q, 4)
@test front(q) == 5
@test dequeue!(q) == 5
@test dequeue!(q) == 6
@test dequeue!(q) == 4
@test isempty(q)
enqueue!(q, 4)
enqueue!(q, 4)
delete!(q, 4)
@test isempty(q)
