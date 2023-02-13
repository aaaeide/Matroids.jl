using IterTools

function rank_and_closure(M, A)
  E, F = M
  if A ⊈ E
    throw(ArgumentError("A is not a subset of E"))
  end

  for (r, Fr) ∈ enumerate(F)
    for closed_set ∈ Fr
      if A ⊆ closed_set
        return (r-1, closed_set) # r-1 since Julia is 1-indexed.
      end
    end
  end
end

function rank(M, A)
  (r, _) = rank_and_closure(M, A)
  return r
end

function closure(M, A)
  (_, cl) = rank_and_closure(M, A)
  return cl
end

function bases(M)
  E, F = M
  r = rank(M, E)

  return family([
    base for base in IterTools.subsets(collect(E), r) 
    if rank(M, base) == r
  ])
end