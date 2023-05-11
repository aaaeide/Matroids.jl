using Graphs, MetaGraphs, StatsBase

include("utils.jl")

"""
    function my_random_erection(n, p, T=UInt16)

Erects a matroid in the same way that Knuth's 1974 matroid erection algorithm does, but only keeps track of the "redundant" sets, ie. the closed sets that are not independent.

n is |E|, p is the list where p[i] is the number of coarsenings (closed sets to increase by 1) at rank i+1. Sets are represented as binary numbers (ø = 0x0000, E = 0xffff when n=16).
"""
function my_random_erection(n, P, T=UInt16, OVERRIDE=[])
  p = copy(P)
  r = 1
  E::T = big"2"^n-1
  
  # R[i] is the set of closed-but-not-independent (redundant) sets of rank i.
  R::Vector{Set{T}} = [Set()]

  # Keeps track of the rank of added closed sets.
  rank = Dict{T, UInt8}(0=>0)

  function generate_covers!()
    for y in R[r]
      t = E-y
      for x in R[r+1]
        # Look for sets in R[r+1] that already contain y and remove excess elements from t.
        if (x&y == y)
          t &= ~x          
          if t == 0 break end
        end
      end
        
      # Insert y ∪ a for all a ∈ t.
      while t > 0
        x = y|(t&-t)
        insert_set!(x)
        t &= ~x
      end
    end
  end

  function insert_set!(x)
    for y in R[r+1] # (+1 due to 1-indexing)
      if haskey(rank, x&y)
        if rank[x&y] < r continue end
        if rank[x&y] == r @goto merge end
        throw(error("got rank > r"))
      end

      if Base.count_ones(x&y) < r continue end

      if find_rank(x&y, r, rank) != false continue end

      @label merge
      # x ∩ y has rank >= r, replace with x ∪ y.
      setdiff!(R[r+1], y)
      insert_set!(x|y)
      return
    end

    push!(R[r+1], x)
    rank[x] = r
  end

  function find_rank(x, r, memo)
    # Checks all sets z in R[1:r]. Returns true if one exists st z ⊆ x.
    if haskey(memo, x) return memo[r] end
    for (i, zi) in enumerate(R[1:r])
      for z in zi
        if x&z == x 
          memo[x] = i-1
          return i-1 
        end
      end
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
      x = reduce(|, [T(1<<(i-1)) for i in sample(1:n, r+1, replace=false)])
      
      for set in R[r]
        if x&set == x continue end
      end
      
      insert_set!(x)
      p[r] -= 1
    end 
    
    # println(r, " ", length(R[r]))
    if E ∈ R[r] break end
    r += 1
  end

  return (n=n, r=r-1, R=R)
end
