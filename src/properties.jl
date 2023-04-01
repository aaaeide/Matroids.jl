# using IterTools

# function rank_and_closure(M, A)
#   E, F = M
#   if A ⊈ E
#     throw(ArgumentError("A is not a subset of E"))
#   end

#   for (r, Fr) ∈ enumerate(F)
#     for closed_set ∈ Fr
#       if A ⊆ closed_set
#         return (r-1, closed_set) # r-1 since Julia is 1-indexed.
#       end
#     end
#   end
# end

# function rank(M, A)
#   (r, _) = rank_and_closure(M, A)
#   return r
# end

# function closure(M, A)
#   (_, cl) = rank_and_closure(M, A)
#   return cl
# end

# function bases(M)
#   E, F = M
#   r = rank(M, E)

#   return family([
#     base for base in IterTools.subsets(collect(E), r) 
#     if rank(M, base) == r
#   ])
# end

# function rank_and_closure(M::KnuthMatroid, A)
#   if A > 2^M.n-1
#     throw(ArgumentError("$(A) is not a subset of $(2^M.n-1)"))
#   end

#   for (r, Fr) in enumerate(M.F)
#     for B in Fr
#       if A & B == A
#         return (r-1, B)
#       end
#     end
#   end
# end

# function bases(M::KnuthMatroid) end