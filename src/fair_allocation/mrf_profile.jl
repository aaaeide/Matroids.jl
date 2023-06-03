using Allocations
import Allocations: bundle_value

"""
    struct MatroidRank <: Profile

A matroid rank valuation profile, representing how each agent values all possible bundles. The profile is constructed from `n` matroids, one for each agents, each matroid over the set of goods [m]. 
"""
struct MatroidRank <: Profile
    matroids::Vector{T} where T <: Matroid
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

"""
    Δ(V::MatroidRank, i::Integer, A::BitMatrix, e::Integer)

Returns the marginal value of adding element e to bundle i in the
allocation A. A is represented as an na x ni BitMatrix, where 
A[i,j] == 1 iff j ∈ A_i.
"""
function Δ(V::MatroidRank, i, A, e)
    Ai´ = copy(A[i, :]); Ai´[e] = 1
    return bv_rank(V.matroids[i], Ai´) - bv_rank(V.matroids[i], A[i, :])
end