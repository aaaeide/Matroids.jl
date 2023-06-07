using Graphs
using Memoize

ground_set(M::ClosedSetsMatroid) = bits_to_set(2^M.n-1)
ground_set(M::FullMatroid) = bits_to_set(2^M.n-1)
ground_set(M::UniformMatroid) = Set(1:M.n)
ground_set(M::GraphicMatroid) = edges(M.g)



"""
    is_indep(M::Matroid, S::Set)

Independence oracle. Determines whether S is independent in M.
"""
function is_indep end

"""
    is_indep(M::ClosedSetsMatroid, S::Integer)

Determines whether a given set S is independent in the matroid M, given by the closed sets of M grouped by rank. Uses (I.1) in Greene (1989).
"""
function is_indep(M::ClosedSetsMatroid, S::Integer)
  t = Base.count_ones(S)

  if t > length(M.F) return false end

  for F in M.F[t]
    if S&F==S return false end
  end

  return true
end

is_indep(M::ClosedSetsMatroid, S) = is_indep(M, set_to_bits(S))


function is_indep(M::FullMatroid, S::Integer)
  card = Base.count_ones(S)
  if card + 1 > length(M.I) return false end
  return S in M.I[card+1]
end

is_indep(M::FullMatroid, S) = is_indep(M, set_to_bits(S))


is_indep(M::UniformMatroid, S) = length(S) <= M.r
is_indep(M::UniformMatroid, S::Integer) = is_indep(M, bits_to_set(S))


function is_indep(m::GraphicMatroid, S)
  edgelist = [e for (i, e) in enumerate(edges(m.g)) if i in S]
  subgraph, _vmap = induced_subgraph(m.g, edgelist)
  return !is_cyclic(subgraph)
end

bv_is_indep(m::Matroid, bv::BitVector) = 
  is_indep(m, [i for (i,e) in enumerate(bv) if e == 1])




"""
    rank(M::Matroid, S)

Returns the rank of the set S in M, ie. the size of the largest independent subset of S.
"""
function rank end


rank(M::Matroid) = M.r


"""
    function rank(M::KnuthMatroid, S::Integer)

Rank ``oracle''. Returns the rank of the set S in the matroid M.
"""
function rank(M::ClosedSetsMatroid, S::Integer)
  for (r, Fr) in enumerate(M.F), B ∈ Fr
      if S&B == S return r-1 end
  end
end
rank(M::ClosedSetsMatroid, S) = rank(M, set_to_bits(S))


function rank(M::FullMatroid, S::Integer)
  for (r, Fr) in enumerate(M.F) for F ∈ Fr
    if S&F == S return r-1 end
  end end
end
rank(M::FullMatroid, S) = rank(M, set_to_bits(S))


rank(M::UniformMatroid, S) = min(length(S), M.r)
rank(M::UniformMatroid, S::Integer) = rank(M, bits_to_set(S))


rank(M::GraphicMatroid, S) = length(S) > 0 ? length(minimal_spanning_subset(M, S)) : 0

"""
    bv_rank(M::Matroid, S::BitVector)

Finds the rank of a subset S of 1:m given as an m-length BitVector bv,
where bv[i] == 1 iff i ∈ S.
"""
bv_rank(M::Matroid, bv::BitVector) = 
  rank(M, [i for (i,e) in enumerate(bv) if e == 1])


"""
    is_circuit(M::Matroid, S)

Determines whether S is a circuit in M (ie. rank(M, S) = |S|-1).
"""
function is_circuit end


"""
    function is_circuit(M::ClosedSetsMatroid, S::Integer)

Determines whether a given set S is a circuit in the matroid M, given by the closed sets of M. Uses (C.1) and (C.2) in Greene (1989).
"""
function is_circuit(M::ClosedSetsMatroid, S::Integer)
  t = Base.count_ones(S)

  for F in M.F[t] # (C.1) S ⊆ F for some F ∈ F_{t-1}.
    if S&F==S @goto C2 end
  end
  return false

  @label C2
  for F in M.F[t-1] # (C.2) |S ∩ F| ≤ r(F) for all F ∈ F_{t-2}.
    if Base.count_ones(S&F) > t-2 return false end
  end

  return true
end


is_circuit(M::FullMatroid, S::Integer) = S in M.C
is_circuit(M::FullMatroid, S) = is_circuit(M, set_to_bits(S))


is_circuit(M::UniformMatroid, S::Integer) = M.r == Base.count_ones(S)-1 <= M.n
is_circuit(M::UniformMatroid, S) = is_circuit(M, set_to_bits(S))


function is_circuit(M::GraphicMatroid, S)
  return rank(M, S) == length(S) - 1
end




"""
    minimal_spanning_subset(M::Matroid, S)

Finds a minimal spanning subset of S in M. If S is the ground set of M, this produces a basis of M.
"""
function minimal_spanning_subset end


"""
    function minimal_spanning_subset(M::ClosedSetsMatroid, A::Integer)

Algorithm 3.1 from Greene (1989). Given a matroid M = (E, F) and some subset A of E, finds a minimal spanning subset of A. If A = E, this finds a basis of M. If A is a basis, this finds A.
"""
minimal_spanning_subset(M::ClosedSetsMatroid, A::Integer) = _mss(M, 0, A)
minimal_spanning_subset(M::ClosedSetsMatroid, A) = minimal_spanning_subset(M, set_to_bits(A))

minimal_spanning_subset(M::FullMatroid, A::Integer) = _mss(M, 0, A)
minimal_spanning_subset(M::FullMatroid, A) = minimal_spanning_subset(M, set_to_bits(A))

function _mss(M, j::Integer, Ā::Integer)
  B = [Ā&F for F in M.F[j+1] if Base.count_ones(Ā&F) > j]

  while length(B) == 0
    if j >= Base.count_ones(Ā)-1 return Ā end
    
    j += 1
    B = [Ā&F for F in M.F[j+1] if Base.count_ones(Ā&F) > j]
  end

  _mss(M, j, Ā&~rand_el(reduce(|, B)))
end


minimal_spanning_subset(M::UniformMatroid, S) = throw("unimplemented")


"""
    minimal_spanning_subset(M::GraphicMatroid, S)

Uses Kruskal's algorithm to find a minimal spanning tree over M.G.
"""
function minimal_spanning_subset(M::GraphicMatroid, S)
  edgelist = [e for (i, e) in enumerate(edges(M.g)) if i in S]
  subgraph, _vmap = induced_subgraph(M.g, edgelist)
  return kruskal_mst(subgraph)
end




"""
    function minimal_spanning_subsets(M::ClosedSetsMatroid, A::Integer)

A modification of Algorithm 3.1 from Greene (1989) that finds all minimal spanning subsets of A ⊆ E, given a matroid M = (E, F). If A = E, this finds the bases of M.
"""
minimal_spanning_subsets(M::ClosedSetsMatroid, A::Integer) = _mss_all(M, 0, A)
minimal_spanning_subsets(M::ClosedSetsMatroid, A) = minimal_spanning_subsets(M, set_to_bits(A))

minimal_spanning_subsets(M::FullMatroid, A::Integer) = _mss_all(M, 0, A)
minimal_spanning_subsets(M::FullMatroid, A) = minimal_spanning_subsets(M, set_to_bits(A))

@memoize function _mss_all(M, j::Integer, Ā::Integer)
  B = [Ā&F for F in M.F[j+1] if Base.count_ones(Ā&F) > j]

  while length(B) == 0
    if j >= Base.count_ones(Ā)-1 return Ā end
    
    j += 1
    B = [Ā&F for F in M.F[j+1] if Base.count_ones(Ā&F) > j]
  end

  bases = Set()
  t = reduce(|, B)
  while t > 0
    x = t&-t
    bases = bases ∪ _mss_all(M, j, Ā&~x)
    t &= ~x
  end
  return bases
end


minimal_spanning_subsets(M::UniformMatroid, A) = throw("unimplemented")
minimal_spanning_subsets(M::GraphicMatroid, A) = throw("unimplemented")




"""
    bases(M::Matroid)

Finds the set of bases of M.
"""
function bases end

bases(M::ClosedSetsMatroid) = _mss_all(M, 0, 2^M.n-1)
bases(M::FullMatroid) = _mss_all(M, 0, 2^M.n-1)


bases(M::UniformMatroid) = throw("unimplemented")
bases(M::GraphicMatroid) = throw("unimplemented")




"""
    closure(M::Matroid, S)

Finds the closure of a set S in M, that, when given a set S ⊆ E, returns the set of elements in x ∈ E such that x can be added to S with no increase in rank. It returns the closed set of the same rank as S, that contains S.
"""
function closure end


function closure(M::ClosedSetsMatroid, S::Integer)
  for Fr in M.F for B ∈ Fr
      if S&B == S return B end
  end end
end
closure(M::ClosedSetsMatroid, S) = closure(M, set_to_bits(S))


function closure(M::FullMatroid, S::Integer)
  for Fr in M.F for B ∈ Fr
      if S&B == S return B end
  end end
end
closure(M::FullMatroid, S) = closure(M, set_to_bits(S))

closure(M::UniformMatroid, S::Integer) = closure(M, bits_to_set(S))
closure(M::UniformMatroid, S) = length(S) < M.r ? S : ground_set(M)


function closure(M::GraphicMatroid, S)
  edgelist = [e for (i, e) in enumerate(edges(M.g)) if i in S]
  _sg, vmap = induced_subgraph(M.g, edgelist)
  return [e for e in edges(M.g) if [e.src, e.dst] ⊆ vmap || e.src == e.dst]
end