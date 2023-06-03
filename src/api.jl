using Graphs
using DataStructures
using Memoize


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
    matroid_partition_knuth_73(Ms, lims=nothing)


Knuth's 1973 Matroid Partitioning algorithm for partitioning a set into subsets independent in various given matroids.

Knuth's description: Given k matroids Ms = [M1, ..., Mk] on the same ground set E, the algorithm finds a k-partitioning [S1, ..., Sk] of the elements of E such that Sj is independent in matroid Mj, and nj <= |Sj| <= nj´, for given limits nj and nj´. 

This implementation drops the upper limit nj´ for each element j (implicitly it is infinity for all matroids). Supply nj in array lims (lims[j] = nj).
"""
function matroid_partition_knuth_73(Ms, lims=nothing)
  n = Ms[1].n; k = length(Ms)
  S0 = Set(1:n) # The unallocated items.
  S = [Set() for _ in 1:k] # The partition-to-be.
  color = Dict(x=>0 for x in 1:n) # color[x] = j iff x ∈ S[j].
  for y in 1:k color[-y] = y end # -y is the 'standard' element of color y.
  succ = [0 for _ in 1:n]

  lims = lims === nothing ? [0 for _ in 1:k] : lims

  function augment(r)
    for x in 1:n succ[x] = 0 end
    
    A = Set(1:n)
    B = r > 0 ? Set(-r) : Set(-j for j in 1:k)
    
    while B != Set()
      C = Set()
      for y ∈ B for x ∈ A
        j = color[y]

        if is_indep(Ms[j], x ∪ setdiff(S[j], y))
          succ[x] = y
          A = setdiff(A, x)
          C = C ∪ x
          if color[x] == 0 repaint(x); return end
        end
      end end
      B = C
    end

    println("$A violates the condition of Theorem 3")
  end

  function repaint(x)
    while x ∈ 1:n
      y = succ[x]
      j = color[x]
      
      if j == 0 setdiff!(S0, x) else setdiff!(S[j], x) end
      
      j = color[y]
      S[j] = S[j] ∪ x
      color[x] = j
      x = y
    end
  end

  # Ensure every part gets at least its lower limit.
  for j in 1:k for _ in 1:lims[j]
    augment(j)
  end end

  # Allocate the rest.
  while S0 != Set()
    augment(0)
  end

  return S
end


# """
#     exchange_graph(Ms::Array{KnuthMatroid}, sets::Array{Integer}, n::Integer)

# Constructs an exchange graph given an m-length array Ms of matroids over the same ground set E, and an m-length array Ss of subsets of E. Each set is represented as an integer of at least n bits, where a 1-bit at bit i denotes that the set contains the element i. Each element is assumed to be in exactly one set. Each set Ss[i] corresponds to the matroid Ms[i], with which the rank can be deduced.

# The exchange graph is a graph whose vertices are the elements of E, and contains the edge (i,j) iff rank_i(S_i - i + j) = rank_i(S_i), where S_i is the set that contains i and rank_i the rank function for the corresponding matroid. In other words, it contains an edge (i,j) iff i can be replaced by j with no loss in rank.

# Returns the exchange graph G and a dictionary d that maps between each element of E and the set in Ss that contains it. 

# NB! Whereas the elements are the set {0, ..., n-1}, the nodes in the graph make up the set {1, ..., n}. That is, vertex v represents element v-1. This is because nodes are 1-indexed in Graphs.jl.
# """
# function exchange_graph(Ms, Ss, n::Integer)
#   G = SimpleDiGraph{UInt8}(n)
#   d = Dict() # e => i, for each e ∈ E, st e ∈ Ss[i].
#   for e in 0:n-1 d[e] = findfirst(S -> 1<<e&S == 1<<e, Ss) end

#   # Checking for each i,j whether element ei can be replaced with ej.
#   for (ei, ej) in Iterators.product(0:n-1, 0:n-1) if ei != ej
#     # Find index of set containing ei.
#     i = d[ei]

#     # Generate ei-ej-replaced set.
#     S´ = Ss[i] & ~(1<<ei) | 1<<ej

#     # Check if rank_i(S_i - ei + ej) = rank_i(S_i)
#     if rank(Ms[i], Ss[i]) == rank(Ms[i], S´)
#       add_edge!(G, ei+1, ej+1) #+1 due to 1-indexing
#     end
#   end end 

#   return (G, d)
# end

"""
Constructs an exchange graph over `m` matroids (E, I_i) and `n` subsets
A_i ⊆ E, given as an n x m BitMatrix A, where A[i,j] = 1 iff element 
j ∈ A_i.

The exchange graph is a graph whose vertices are the elements of E, and contains the edge (i,j) iff rank_i(S_i - i + j) = rank_i(S_i), where S_i is the set that contains i and rank_i the rank function for the corresponding matroid. In other words, it contains an edge (i,j) iff i can be replaced by j with no loss in rank.
"""
function exchange_graph(Ms::Vector{T}, A::BitMatrix) where T <: Matroid
  n, m = size(A)
  @assert n == length(Ms); @assert m == Ms[1].n; @assert sum(A) == m
  D = SimpleDiGraph(m)
  
  # Checking for each i,j whether element ei can be replaced with ej.
  for ei in 1:m, ej in setdiff(1:m, ei) 
    # Find index of set containing ei.
    i = findfirst(==(1), A[:, ei])

    # Generate ei-ej-replaced set.
    Ai´ = A[i, :]
    Ai´[ei] = 0; Ai´[ej] = 1

    # Check if rank_i(A_i - ei + ej) = rank_i(A_i)
    if bv_rank(Ms[i], Ai´) == bv_rank(Ms[i], A[i, :])
      add_edge!(D, ei, ej)
    end
  end

  return D
end

"""
    find_shortest_path(D, from, to)

Finds a shortest path [s,...,t], where s ∈ from and t ∈ to.
"""
function find_shortest_path(D, from, to)
  X = intersect(from, to)
  if length(X) > 0
    return [X[1]]
  end

  ds = dijkstra_shortest_paths(D, from)
  paths = []

  for g in to
    path = [ds.parents[g], g]
    
    while path[1] ∉ from
      if path[1] == 0 @goto skip end
      pushfirst!(path, ds.parents[path[1]])
    end

    push!(paths, path)
    @label skip
  end

  if length(paths) == 0 return nothing end
  return argmin(length, paths)
end


"""
    transfer!(A, i, path, D)

Creates a new allocation from A, where the goods have between
transferred along path, path[1] ending up in agent i's bundle. 

Returns the augmented allocation and an updated exchange graph D.
"""
function transfer!(A, i, path)
  # At every iteration, x receives the next good in the path.
  x = i
  for g in path
    # y is the current owner, who loses g.
    y = findfirst(==(1), A[:, g])
    A[x, g] = 1; A[y, g] = 0
    x = y
  end

  @assert sum(A) == size(A)[2]
  return A
end


# """
#     function find_transfer_path(G, d, Ms, Ss, i::Integer, j::Integer, n::Integer)

# Accepts an m-length array Ms of matroids over a ground set E, an m-length array Ss of subsets of E, the size of the universe |E| = n, and indexes 0 <= i <= m, 0 <= j <= m. Finds the shortest **transfer path** between i and j.

# A transfer path is a path (s, e_1, ..., e_k) in the exchange graph G, such that 
# - Δ(Ms[i], Ss[i], e_1) = 1,
# - rank(Ms[x], Ss[x] - e_x + e_{x+1}) = rank(Ms[x], Ss[x]) for each x in 2:k-1.
# """
# function find_transfer_path(G, d, Ms, Ss, n::Integer, i::Integer, j::Integer)
#   # Compute the set F_i of elements with positive marginal gain for i.
#   # TODO: Remember this between iterations.
#   F_i = [e for e in vertices(G) if Δ(Ms[i], Ss[i], e-1) == 1] 
  
#   add_vertex!(G)
#   s = last(vertices(G))
#   for e in F_i add_edge!(G, s, e) end

#   # This BFS gives us the path [s, e_1, ..., e_k] of elements (0..n-1) or false.
#   return bfs(G, s, v -> d[v-1] == j)
# end


# function transfer!(i, p, Ms, Ss, d, n)
#   @assert p[1] == n+1 # The source node, not an element in E.

#   # Ss[i] receives p[2] for a marginal gain of 1.
#   Ss[i] |= 1<<(p[2]-1)
  
#   # For each path position j in 2:end-1, d[j] should lose p[j] and get p[j+1].
#   for j in 2:length(p)-1
#     r = rank(Ms[d[j]], Ss[d[j]])
#     Ss[d[j]] &= ~1<<(p[j]-1) # Lose an element.
#     Ss[d[j]] |= 1<<(p[j+1]-1) # Gain next element.
#     r´ = rank(Ms[d[j]], Ss[d[j]])
#     @assert r == r´
#   end

#   # The last set simply loses an element.
#   Ss[d[p[end]-1]] &= ~1<<(p[end]-1)
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