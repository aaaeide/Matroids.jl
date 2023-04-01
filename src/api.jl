using Graphs
using DataStructures

"""
    rank(M::KnuthMatroid, S::Integer)

Rank ``oracle''. Returns the rank of the set S in the matroid M.
"""
function rank(M::KnuthMatroid, S::Integer)
  for (r, Fr) in enumerate(M.F)
    for B ∈ Fr
      if S&B == S return r-1 end
    end
  end
end

function rank(M::KnuthMatroid, S::Set) return rank(M, set_to_bits(S)) end


"""
    exchange_graph(Ms::Array{KnuthMatroid}, sets::Array{Integer}, n::Integer)

Constructs an exchange graph given an m-length array Ms of matroids over the same ground set E, and an m-length array Ss of subsets of E (represented as integers of at least n bits). Each element is assumed to be in exactly one set. Each set Ss[i] corresponds to the matroid Ms[i], with which the rank can be deduced.

The exchange graph is a graph whose vertices are the elements of E, and contains the edge (i,j) iff rank_i(S_i - i + j) = rank_i(S_i), where S_i is the set that contains i and rank_i the rank function for the corresponding matroid. In other words, it contains an edge (i,j) if i can be replaced by j with no loss in rank.
"""
function exchange_graph(Ms, Ss, n::Integer)
  G = SimpleDiGraph{UInt8}(n)

  # TODO: Consider using MetaGraphs to store which element a vertex represents.

  # Checking for each i,j whether element ei can be replaced with ej.
  for (ei, ej) in Iterators.product(0:n-1, 0:n-1) if ei != ej
    # Find index of set containing ei.
    i = findfirst(S -> 1<<ei&S == 1<<ei, Ss)

    # Generate ei-ej-replaced set.
    S´ = Ss[i] & ~(1<<ei) | 1<<ej

    # Check if rank_i(S_i - ei + ej) = rank_i(S_i)
    if rank(Ms[i], Ss[i]) == rank(Ms[i], S´)
      add_edge!(G, ei+1, ej+1) #+1 due to 1-indexing
    end
  end end 

  return G
end

"""
    function Δ(M, S, e)

Returns the marginal gain of adding e to the set S, given the matroid M.

0 ≤ e ≤ M.n-1
"""
function Δ(M, S, e)
  return rank(M, S | 1<<e) - rank(M, S)
end

function find_transfer_path(Ms, Ss, i::Integer, j::Integer, n::Integer)
  G = exchange_graph(Ms, Ss, n)
  
  # Compute the set F_i of elements with positive marginal gain for i.
  # TODO: Remember this between iterations.
  F_i = [e for e in vertices(G) if Δ(Ms[i], Ss[i], e-1) == 1]
  
  # Add source node.
  add_vertex!(G)
  s = last(vertices(G))
  @assert neighbors(G, s) == []

  # Add edge (s,e) for each e ∈ F_i.
  for e in F_i @assert add_edge!(G, s, e) end

  # Find shortest path from s to j
  q = Queue{Integer}()
  parent = Dict()
  enqueue!(q, s)
  while length(q) > 0
    v = dequeue!(q)
    if 1<<(v-1)&Ss[j] == 1<<(v-1) break end # TODO
    for u in neighbors(G, v) 
      parent[u] = v
      enqueue!(q, u)
    end
  end

  x = y
end