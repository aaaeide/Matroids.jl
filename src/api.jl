using Graphs
using DataStructures
# using Memoize

include("utils.jl")

"""
    function is_indep(M::ClosedSetsMatroid, S::Integer)

Determines whether a given set S is independent in the matroid M, given by the closed sets of M grouped by rank. Uses (I.1) in Greene (1989).
"""
function is_indep(M::ClosedSetsMatroid, S::Integer)
  t = Base.count_ones(S)

  for F in M.F[t]
    if Base.count_ones(S&F) > t-1 return false end
  end

  return true
end

function is_indep(M::FullMatroid, S::Integer)
  card = Base.count_ones(S)
  return S in M.I[card+1]
end

function is_indep(M::UniformMatroid, S::Integer)
  return Base.count_ones(S) <= M.r
end

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

"""
    function minimal_spanning_subsets(M::ClosedSetsMatroid, A::Integer)

Algorithm 3.1 from Greene (1989). Given a matroid M = (E, F) and some subset A of E, finds a minimal spanning subset of A. If A = E, this finds a basis of M.
"""
minimal_spanning_subset(M::ClosedSetsMatroid, A::Integer) = _mss(M, 0, A)

function _mss(M::ClosedSetsMatroid, j::Integer, Ā::Integer)
  B = [Ā&F for F in M.F[j+1] if Base.count_ones(Ā&F) > j]

  while length(B) == 0
    if j >= Base.count_ones(Ā)-1 return Ā end
    
    j += 1
    B = [Ā&F for F in M.F[j+1] if Base.count_ones(Ā&F) > j]
  end

  _mss(M, j, Ā&~rand_el(reduce(|, B)))
end

"""
    function minimal_spanning_subsets(M::ClosedSetsMatroid, A::Integer)

A modification of Algorithm 3.1 from Greene (1989) that finds all minimal spanning subsets of A ⊆ E, given a matroid M = (E, F). If A = E, this finds the bases of M.
"""
minimal_spanning_subsets(M::ClosedSetsMatroid, A::Integer) = _mss_all(M, 0, A)

# TODO: @memoize
function _mss_all(M::ClosedSetsMatroid, j::Integer, Ā::Integer)
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

bases(M::ClosedSetsMatroid) = _mss_all(M, 0, 2^M.n-1)

"""
    function rank(M::KnuthMatroid, S::Integer)

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

function rank(M::UniformMatroid, S::Integer)
  return min(Base.count_ones(S), M.r)
end
  
function rank(M, S) return rank(M, set_to_bits(S)) end


"""
    function Δ(M, S, e)

Returns the marginal gain of adding e to the set S, given the matroid M.

0 ≤ e ≤ M.n-1
"""
function Δ(M, S, e)
  return rank(M, S | 1<<e) - rank(M, S)
end


"""
    function matroid_partition(Ms::Vector{Matroid})

Partitions a set of elements into k subsets independent in the k given matroids Ms. The matroids are assumed to be over the same ground set E.

Wants to create an "even" partition; will at each step enlarge the set Ss[i], where Ss[i] is not a basis and is of smallest cardinality among the sets in Ss.
"""
function matroid_partition(Ms::Vector{Matroid})
  k = length(Ms); n = Ms[1].n

  # Initially, all k subsets are empty, and the "pot" is full.
  Ss = push!(zeros(UInt16, k), 2^n-1)
  push!(Ms, FreeMatroid(n))
  
  # can_grow[i] == true iff Ss[i] is not a basis in Ms[i].
  can_grow = trues(k)
  while true in can_grow
    # Find the index i of the first set of least cardinality.
    smallest = n+1; i = -1
    for p in 1:k if can_grow[p]
      card = Base.count_ones(Ss[p])
      if card < smallest
        smallest = card; i = p
      end
    end end

    Ss |> pb
    readline()

    G, d = exchange_graph(Ms, Ss, n)
    path = find_transfer_path(G, d, Ms, Ss, n, i, k+1)
    if path == false
       can_grow[i] = false
    else
      transfer!(i, path, Ms, Ss, d, n)
    end
  end

  popat!(Ms, k+1) # Cleanup.
  return Ss
end


"""
    exchange_graph(Ms::Array{KnuthMatroid}, sets::Array{Integer}, n::Integer)

Constructs an exchange graph given an m-length array Ms of matroids over the same ground set E, and an m-length array Ss of subsets of E. Each set is represented as an integer of at least n bits, where a 1-bit at bit i denotes that the set contains the element i. Each element is assumed to be in exactly one set. Each set Ss[i] corresponds to the matroid Ms[i], with which the rank can be deduced.

The exchange graph is a graph whose vertices are the elements of E, and contains the edge (i,j) iff rank_i(S_i - i + j) = rank_i(S_i), where S_i is the set that contains i and rank_i the rank function for the corresponding matroid. In other words, it contains an edge (i,j) iff i can be replaced by j with no loss in rank.

Returns the exchange graph G and a dictionary d that maps between each element of E and the set in Ss that contains it. 

NB! Whereas the elements are the set {0, ..., n-1}, the nodes in the graph make up the set {1, ..., n}. That is, vertex v represents element v-1. This is because nodes are 1-indexed in Graphs.jl.
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
    function find_transfer_path(G, d, Ms, Ss, i::Integer, j::Integer, n::Integer)

Accepts an m-length array Ms of matroids over a ground set E, an m-length array Ss of subsets of E, the size of the universe |E| = n, and indexes 0 <= i <= m, 0 <= j <= m. Finds the shortest **transfer path** between i and j.

A transfer path is a path (s, e_1, ..., e_k) in the exchange graph G, such that 
- Δ(Ms[i], Ss[i], e_1) = 1,
- rank(Ms[x], Ss[x] - e_x + e_{x+1}) = rank(Ms[x], Ss[x]) for each x in 2:k-1.
"""
function find_transfer_path(G, d, Ms, Ss, n::Integer, i::Integer, j::Integer)
  # Compute the set F_i of elements with positive marginal gain for i.
  # TODO: Remember this between iterations.
  F_i = [e for e in vertices(G) if Δ(Ms[i], Ss[i], e-1) == 1] 
  
  add_vertex!(G)
  s = last(vertices(G))
  for e in F_i add_edge!(G, s, e) end

  # This BFS gives us the path [s, e_1, ..., e_k] of elements (0..n-1) or false.
  return bfs(G, s, v -> d[v-1] == j)
end


function transfer!(i, p, Ms, Ss, d, n)
  @assert p[1] == n+1 # The source node, not an element in E.

  # Ss[i] receives p[2] for a marginal gain of 1.
  Ss[i] |= 1<<(p[2]-1)
  
  # For each path position j in 2:end-1, d[j] should lose p[j] and get p[j+1].
  for j in 2:length(p)-1
    r = rank(Ms[d[j]], Ss[d[j]])
    Ss[d[j]] &= ~1<<(p[j]-1) # Lose an element.
    Ss[d[j]] |= 1<<(p[j+1]-1) # Gain next element.
    r´ = rank(Ms[d[j]], Ss[d[j]])
    @assert r == r´
  end

  # The last set simply loses an element.
  Ss[d[p[end]-1]] &= ~1<<(p[end]-1)
end

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