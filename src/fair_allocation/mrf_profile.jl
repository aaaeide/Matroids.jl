using Allocations

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
value(V::MatroidRank, i, S::BitVector) = bv_rank(V.matroids[i], S)
value(V::MatroidRank, i, g::Int) = value(V, i, Set(g))

# Disambiguation:
"""
    value(V::MatroidRank, i, A::Allocation; indep=false)

Returns the agent i's realized value in allocation A. Supply 
indep=true if you know A_i is independent in M_i, to save a call
to the rank function.
"""
value(V::MatroidRank, i, A::Allocation; indep=false) = 
    indep ? length(bundle(A, i)) : value(V, i, bundle(A, i))

"""
    Δ(V::MatroidRank, A, i, g)

Returns the marginal value of adding element g to bundle A_i.
"""
Δ(V::MatroidRank, A, i, g) = 
    value(V, i, bundle(A, i) ∪ g) - value(V, i, bundle(A, i))

"""
    Δ(V::MatroidRank, i::Integer, A::BitMatrix, g::Integer)

Returns the marginal value of adding element g to bundle A_i. A is represented as an na x ni BitMatrix, where 
A[i,j] == 1 iff j ∈ A_i.
"""
function Δ(V::MatroidRank, i::Integer, A::BitMatrix, g::Integer)
    Ai´ = copy(A[i, :]); Ai´[g] = 1
    return bv_rank(V.matroids[i], Ai´) - bv_rank(V.matroids[i], A[i, :])
end