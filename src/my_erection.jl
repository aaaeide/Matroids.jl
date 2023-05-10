using Graphs, MetaGraphs, StatsBase

include("utils.jl")

function VERBOS(thresh) thresh >= 100 end
function bst(x) bitstring(x)[1:end] end

"""
    function my_random_erection(n, p, T=UInt16)

Erects a matroid in the same way that Knuth's 1974 matroid erection algorithm does, but only keeps track of the "redundant" sets, ie. the closed sets that are not independent.

n is |E|, p is the list where p[i] is the number of coarsenings (closed sets to increase by 1) at rank i+1. Sets are represented as binary numbers (ø = 0x0000, E = 0xffff when n=16).
"""
function my_random_erection(n, P, T=UInt16, OVERRIDE=[])
  p = copy(P)
  r = 1
  E::T = big"2"^n-1

  log = []
  
  # R[i] is the set of closed-but-not-independent (redundant) sets of rank i.
  R::Vector{Set{T}} = [Set()]

  # Keeps track of the rank of added closed sets.
  rank = Dict{T, UInt8}(0=>0)

  # # G contains the edge (X, Y) if X is a cover ("descendant") of Y. Then, X ⊇ Y.
  # G = MetaDiGraph()

  # # Each node is labeled with the bitstring of the set it represents.
  # set_indexing_prop!(G, :label) 

  # function G_get(x) 
  #   try
  #     G[bitstring(x), :label] 
  #   catch
  #     return nothing
  #   end
  # end

  # function G_set(x)
  #   if G_get(x) === nothing
  #     add_vertex!(G)
  #     set_prop!(G, nv(G), :label, bitstring(x))
  #     return true
  #   end

  #   return false
  # end

  # function G_neighbs(x)
  #   if G_get(x) === nothing
  #     return []
  #   end

  #   return neighbors(G, G_get(x))
  # end

  # function G_neighbs(xs...)
  #   return reduce(vcat, [G_neighbs(x) for x in xs])
  # end


  function generate_covers!()
    for y in R[r]
      VERBOS(2) && println("\n=== generating covers for $(bst(y)) ===")
      t = E-y
      for x in R[r+1]
        VERBOS(0) && print("\n\tcomparing with \t$(bst(x)) ")
        # Look for sets in R[r+1] that already contain y.
        if (x&y == y)
          VERBOS(0) && print("*")
          # Remove excess elements from t.
          t &= ~x

          # Add edge (x, y) to remember that x ⊇ y.
          # add_edge!(G, (G_get(x), G_get(y)))
          
          if t == 0 break end
        end
      end

      VERBOS(1) && println("\n\tadding all \t$(bst(t))\n")
        
      # Insert y ∪ a for all a ∈ t.
      while t > 0
        x = y|(t&-t)
        insert_set!(x)#, [G_get(y)])
        t &= ~x
      end
    end
  end

  function insert_set!(x)#, parents)
    VERBOS(1) && readline()
    VERBOS(2) && println(" trying to add $(bst(x)) to rank $r")

    for y in R[r+1] # (+1 due to 1-indexing)
      if haskey(rank, x&y)
        if rank[x&y] < r 
          VERBOS(0) && println("\ttable ok ($(rank[x&y])): $(bst(y))")
          continue 
        end
        if rank[x&y] == r @goto merge end
        throw(error("got rank > r"))
      end

      if Base.count_ones(x&y) < r
        VERBOS(0) && println("\tcardinality ok: $(bst(y))")
        continue
      end

      # parents = vcat(parents, G_neighbs(y))
      # if check_rank(x&y, parents)
      #   VERBOS(0) && println("RANK CHECK FIRE")
      #   continue
      # end

      if check_rank(x&y, r) continue end

      @label merge
      VERBOS(2) && println("\t$(bst(x)) ∩ $(bst(y)) = $(bst(x&y)) has rank == $r")
      VERBOS(2) && println("\treplacing with $(bst(x|y))")

      # x ∩ y has rank >= r, replace with x ∪ y.
      setdiff!(R[r+1], y)
      # parents = vcat(parents, G_neighbs(y))
      # rem_vertex!(G, G_get(x)); rem_vertex!(G, G_get(y))
      insert_set!(x|y)#, parents)
      return
    end

    VERBOS(3) && println("\tadding $(bst(x)) to rank $r")

    push!(R[r+1], x)
    rank[x] = r
    # G_set(x)
    # for p in parents add_edge!(G, (G_get(x), p)) end
  end

  # function check_rank(x, fringe)
  #   # Runs a BFS, starting with the nodes of G in fringe, looking
  #   # for a set y st. x ⊆ y. if x is contained in a closed set (of lower rank) return true, else false.
  #   VERBOS(1) && print("\tgraph check: $(bst(x))")
  #   visited = Set()
  #   while length(fringe) > 0
  #     y = parse(T, get_prop(G, popfirst!(fringe), :label), base=2)
  #     if y in visited continue end

  #     VERBOS(0) && print("\n\t\t     $(bst(y))")

  #     if x&y == x 
  #       VERBOS(1) && print(" - ok ($(rank[y]))\n")
  #       return true
  #     end
      
  #     fringe = vcat(fringe, G_neighbs(y))
  #     push!(visited, y)
  #   end

  #   VERBOS(1) && print(" - not ok, merge\n")
  #   return false
  # end

  function check_rank(x, r)
    # Checks all sets z in R[r]. Returns true if one exists st z ⊆ x.
    for z in R[r]
      if x&z == x 
        VERBOS(0) && println("\tmanual check ok: $(bst(y))")
        return true 
      end
    end
    return false
  end
  

  while r < n
    # Redundant sets of rank r+1 will go here.
    push!(R, Set())
    generate_covers!()

    if length(OVERRIDE) > 0
      if r < length(OVERRIDE)
        for set in OVERRIDE[r]
          insert_set!(set)#, [])
        end
      end
    else
      # Apply coarsening.
      push!(log, [])
      while r <= length(p) && p[r] > 0 && E ∉ R[r+1]
        # Generate a random set of cardinality r+1.
        VERBOS(2) && println("\n********      coarsening      ********")
        x = reduce(|, [T(1<<(i-1)) for i in sample(1:n, r+1, replace=false)])
        
        for set in R[r]
          if x&set == x continue end
        end

        push!(log[r], x)
        
        insert_set!(x)#, [])
        p[r] -= 1
      end 
    end
      
    VERBOS(2) && println("\nnext rank.")#\nR[$(r+1)] = $(R[r+1])\n")

    if E ∈ R[r] break end
    r += 1
  end

  # println(log)
  return (n=n, r=r-1, R=R)#, G=G)
end
