using DataStructures

head = HalfPath()
tail = HalfPath()
for x in 3:-1:1; enqueue!(head, x); end
for x in 4:6; enqueue!(tail, x); end
path = Path(head, tail)
@test collect(path) == [x for x in 1:6]
