KnuthMatroid = Tuple(UInt16, Vector{Set{UInt16}})

function diff(A, B)
  A & ~B
end

function subseteq(A, B)
  A & B == A
end

