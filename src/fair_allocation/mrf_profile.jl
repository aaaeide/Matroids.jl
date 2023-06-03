using Allocations
import Allocations: bundle_value

"""
    struct MatroidRank <: Profile

A matroid rank valuation profile, representing how each agent values all possible bundles. The profile is constructed from `n` matroids, one for each agents, each matroid over the set of goods [m]. 
"""
struct MatroidRank <: Profile
    matroids::Vector{Matroid}
    m::Int
end

na(V::MatroidRank) = length(V.matroids)
ni(V::MatroidRank) = V.m

value(V::MatroidRank, i, S) = rank(V.matroids[i], S)
value(V::MatroidRank, i, S::BitVector) = 
    rank(V.matroids[i], [g for g in 1:ni(V) if S[g] == 1])
value(V::MatroidRank, i, g::Int) = value(V, i, Set(g))

# Disambiguation:
value(V::MatroidRank, i, A::Allocation) = bundle_value(V, i, A)

Δ(V::MatroidRank, i, A, e) = Δ(V.matroids[i], A, e)