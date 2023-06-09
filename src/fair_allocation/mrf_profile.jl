using Allocations

"""
    struct MatroidRank <: Profile

A matroid rank valuation profile, representing how each agent values all possible bundles. The profile is constructed from `n` matroids, one for each agents, each matroid over the set of goods [m]. 
"""
struct MatroidRank <: Profile
    Ms::Vector{T} where T <: Matroid
    m::Int
end

Allocations.na(V::MatroidRank) = length(V.Ms)
Allocations.ni(V::MatroidRank) = V.m

Allocations.value(V::MatroidRank, i, S) = rank(V.Ms[i], S)
Allocations.value(V::MatroidRank, i, g::Int) = value(V, i, Set(g))

Allocations.value(V::MatroidRank, i, A::Allocation) = value(V, i, bundle(A, i))

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
    return bv_rank(V.Ms[i], Ai´) - bv_rank(V.Ms[i], A[i, :])
end


"""
Generate a MatroidRank valuation profile with n matroids over m elements. 
"""
function gen_matroidrank_profile(n, gen_matroid, T)
  function gen()
    ms = Array{T}(undef, n)

    Threads.@threads for i in 1:n
      ms[i] = gen_matroid()
    end

    return MatroidRank(ms, ms[1].n)
  end

  return gen
end

function random_er_graph(m)
  n = rand(trunc(Integer, m/2):3m)
  return erdos_renyi(n, m)
end

function random_ws_graph(m)
  n = rand([x for x in 2:ceil(Integer, sqrt(2m)) if 2m%x == 0 && iseven(x)])
  k = 2m ÷ n
  (k, n) = sort([n,k])

  return watts_strogatz(n, k, rand())
end

function random_ba_graph(m)
  k = rand([x for x in 1:ceil(Integer, sqrt(m)) if m%x == 0])
  n = m ÷ k + k

  return barabasi_albert(n, k)
end
