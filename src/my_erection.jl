using Graphs, MetaGraphs, StatsBase

include("utils.jl")

function VERBOSE(thresh) thresh >= 100 end
function bst(x) bitstring(x)[1:end] end

"""
    function my_random_erection(n, p, T=UInt16)

Erects a matroid in the same way that Knuth's 1974 matroid erection algorithm does, but only keeps track of the "redundant" sets, ie. the closed sets that are not independent.

n is |E|, p is the list where p[i] is the number of coarsenings (closed sets to increase by 1) at rank i+1. Sets are represented as binary numbers (ø = 0x0000, E = 0xffff when n=16).
"""
function my_random_erection(n, p, T=UInt16)
  r = 1
  E::T = big"2"^n-1
  
  # R[i] is the set of closed-but-not-independent (redundant) sets of rank i.
  R::Vector{Set{T}} = [Set()]

  # Keeps track of the rank of added closed sets.
  rank = Dict{T, UInt8}()

  # G contains the edge (X, Y) if X is a cover ("descendant") of Y. Then, X ⊇ Y.
  G = MetaDiGraph()

  # Each node is labeled with the bitstring of the set it represents.
  set_indexing_prop!(G, :label) 

  function G_get(x) 
    try
      G[bitstring(x), :label] 
    catch
      return nothing
    end
  end

  function G_set(x)
    if G_get(x) === nothing
      add_vertex!(G)
      set_prop!(G, nv(G), :label, bitstring(x))
      return true
    end

    return false
  end

  function G_neighbs(x)
    if G_get(x) === nothing
      return []
    end

    return neighbors(G, G_get(x))
  end

  function G_neighbs(xs...)
    return reduce(vcat, [G_neighbs(x) for x in xs])
  end


  function generate_covers!()
    for y in R[r]
      VERBOSE(2) && println("\n=== generating covers for $(bst(y)) ===")
      t = E-y
      for x in R[r+1]
        VERBOSE(1) && print("\n\tcomparing with \t$(bst(x)) ")
        # Look for sets in R[r+1] that already contain y.
        if (x&y == y)
          VERBOSE(1) && print("*")
          # Remove excess elements from t.
          t &= ~x

          # Add edge (x, y) to remember that x ⊇ y.
          add_edge!(G, (G_get(x), G_get(y)))
          
          if t == 0 break end
        end
      end

      VERBOSE(1) && println("\n\tadding all \t$(bst(t))")
        
      # Insert y ∪ a for all a ∈ t.
      while t > 0
        x = y|(t&-t)
        insert_set!(x, [G_get(y)])
        t &= ~x
      end
    end
  end

  function insert_set!(x, parents)
    VERBOSE(1) && readline()
    VERBOSE(2) && print(" trying to add $(bst(x)) to rank $r")
    VERBOSE(2) && if Base.count_ones(x) > r print(" (REDUNDANCE)\n") else print("\n") end

    # We are trying to insert x as a closed set of rank r.
    for y in R[r+1] # (+1 due to 1-indexing)
      if haskey(rank, x&y)
        if rank[x&y] < r continue end
        if rank[x&y] == r @goto merge end
        throw(error("got rank > r"))
      end

      # TODO: Crashes here
      parents = G_neighbs(x, y)
      if check_rank(x&y, parents) < r
        continue
      end

      @label merge
      VERBOSE(2) && println("\t$(bst(x)) ∩ $(bst(y)) = $(bst(x&y)) has rank == $r")
      VERBOSE(2) && println("\treplacing with $(bst(x|y))")

      # x ∩ y has rank >= r, replace with x ∪ y.
      setdiff!(R[r+1], y)
      parents = vcat(parents, G_neighbs(y))
      # rem_vertex!(G, G_get(x)); rem_vertex!(G, G_get(y))
      insert_set!(x|y, parents)
    end

    VERBOSE(3) && println("\tadding $(bst(x)) to rank $r")

    push!(R[r+1], x)
    rank[x] = r
    G_set(x)
    for p in parents add_edge!(G, (G_get(x), p)) end
  end

  function check_rank(x, fringe)
    # Runs a BFS, starting with the nodes of G in fringe, looking
    # for a set y st. x ⊆ y. Returns the rank of y, or false.
    visited = Set()
    while length(fringe) > 0
      y = parse(T, get_prop(G, popfirst!(fringe), :label), base=2)
      if y in visited continue end

      if x&y == x return rank[y] end
      
      fringe = vcat(fringe, G_neighbs(y))
      push!(visited, y)
    end

    return false
  end
  

  while r < n
    # Redundant sets of rank r+1 will go here.
    push!(R, Set())
    generate_covers!()

    # Apply coarsening.
    while r <= length(p) && p[r] > 0 && E ∉ R[r+1]
      # Generate a random set of cardinality r+1.
      VERBOSE(2) && println("\n********      coarsening      ********")
      x = reduce(|, [T(1<<(i-1)) for i in sample(1:n, r+1, replace=false)])
      insert_set!(x, [])
      p[r] -= 1
    end

    # if E ∈ R[r] break end
    r += 1
  end

  return (n=n, r=r-1, R=R)
end
