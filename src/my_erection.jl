include("utils.jl")

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
  rank = Dict{T, UInt8}()
  
  while r < n
    # Redundant sets of rank r+1 will go here.
    push!(R, Set())

    # Add covers of r-rank redundant sets.
    for y in R[r]
      t = E - y
      
      for x in R[r+1]
        if (x&y == y) t&= ~x end
      end
      
      while t > 0
        x = y|(t&-t)
        println(bitstring(x))
        t &= ~x
      end
    end

    if r <= length(p)
      @assert Base.count_ones(s) > r "Rank cannot exceed cardinality!"
      
      push!(R[r+1], s) # rank table add all subsets of rank r+1
    end

    # if E ∈ R[r] break end
    r += 1
    readline()
  end
end

"""
    function push_set!(s, R, r, rank)

Adds set s to R[r+1], augmenting s if necessary to ensure that no two sets in R[r+1] have an intersection of rank greater than r.
"""
function push_set!(s, R, r, rank)
  for t in R[r+1]
    
  end
end