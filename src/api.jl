using Graphs
using DataStructures

include("types.jl")

"""
    rank(M::KnuthMatroid, S::Integer)

Rank ``oracle''. Returns the rank of the set S in the matroid M.
"""
function rank(M::ClosedSetsMatroid, S::Integer)
  for (r, Fr) in enumerate(M.F)
    for B ∈ Fr
      if S&B == S return r-1 end
    end
  end
end


function rank(M::FullMatroid, S::Integer)
  for (r, Ir) in enumerate(M.I)
    if S in Ir return r-1 end
  end
end

function rank(M::FreeMatroid, S::Integer)
  return Base.count_ones(S)
end
  
function rank(M, S::Set) return rank(M, set_to_bits(S)) end

"""
    exchange_graph(Ms::Array{KnuthMatroid}, sets::Array{Integer}, n::Integer)

Constructs an exchange graph given an m-length array Ms of matroids over the same ground set E, and an m-length array Ss of subsets of E (represented as integers of at least n bits). Each element is assumed to be in exactly one set. Each set Ss[i] corresponds to the matroid Ms[i], with which the rank can be deduced.

The exchange graph is a graph whose vertices are the elements of E, and contains the edge (i,j) iff rank_i(S_i - i + j) = rank_i(S_i), where S_i is the set that contains i and rank_i the rank function for the corresponding matroid. In other words, it contains an edge (i,j) if i can be replaced by j with no loss in rank.

Returns the exchange graph G and a dictionary d that maps between each element of E and the set in Ss that contains it.
"""
function exchange_graph(Ms, Ss, n::Integer)
  G = SimpleDiGraph{UInt8}(n)
  d = Dict() # e => i, for each e ∈ E, st e ∈ Ss[i].
  for e in 0:n-1 d[e] = findfirst(S -> 1<<e&S == 1<<e, Ss) end

  # Checking for each i,j whether element ei can be replaced with ej.
  for (ei, ej) in Iterators.product(0:n-1, 0:n-1) if ei != ej
    # Find index of set containing ei.
    i = d[ei]

    # Generate ei-ej-replaced set.
    S´ = Ss[i] & ~(1<<ei) | 1<<ej

    # Check if rank_i(S_i - ei + ej) = rank_i(S_i)
    if rank(Ms[i], Ss[i]) == rank(Ms[i], S´)
      add_edge!(G, ei+1, ej+1) #+1 due to 1-indexing
    end
  end end 

  return (G, d)
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
  G, d = exchange_graph(Ms, Ss, n)

  
  # Compute the set F_i of elements with positive marginal gain for i.
  # TODO: Remember this between iterations.
  F_i = [e for e in vertices(G) if Δ(Ms[i], Ss[i], e-1) == 1]
  println("F_i: $F_i")
  
  # Add source node.
  add_vertex!(G)
  s = last(vertices(G))
  @assert neighbors(G, s) == []

  
  # Add edge (s,e) for each e ∈ F_i.
  for e in F_i @assert add_edge!(G, s, e) end
  println(neighbors(G, s))
  
  nodelabel = vcat([x for x in 0:n-1], "s")
  draw_graph(G, nodelabel)

  println("FINDING SHORTEST PATH BETWEEN $s AND $(bitstring(Ss[j]))")

  return G, bfs(G, s, v -> d[v] == j)
end


# function transfer(i, path, G, d, n)
#   @assert path[1] == n # The source node, not an element.
#   for e in path[2:end]
    
#   end
# end

"""
    function bfs(G, s, fn)
  
Performs a BFS over the graph G starting from s, until it finds the first vertex t such that fn(t) === true. Then returns the path (s,...,t).

Returns false if no such vertex t exists.
"""
function bfs(G, s, fn)
  Q = Queue{Integer}()
  parent = Dict()
  seen = Set()
  t = nothing

  for v in neighbors(G, s) 
    enqueue!(Q, v)
    parent[v] = s 
  end

  while length(Q) > 0
    v = dequeue!(Q)
    if v in seen continue end
    if fn(v) t = v; break end
    for u in neighbors(G, v)
      enqueue!(Q, u)
      if !haskey(parent, u) parent[u] = v end
    end
    push!(seen, v)
  end

  if t === nothing return false end

  path = [t]
  while parent[t] !== s
    push!(path, parent[t])
    t = parent[t]
  end
  push!(path, s)
  return reverse(path)
end