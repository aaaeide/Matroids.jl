include("erect.jl")
include("types.jl")

function randomized_erect_v1(n, p, T=UInt16)
  # Initialize.
  r = 1
  pr = 0
  F::Vector{Set{T}} = [Set(T(0))]
  E::T = big"2"^n-1
  rank = Dict{T, UInt8}()

  # Populate rank table with 100+cardinality for all subsets of E.
  k=1; rank[0]=100;
  while (k<=E)
    for i in 0:k-1 rank[k+i] = rank[i]+1 end
    k=k+k;
  end

  F = [Set(0)] # F[r] is the family of closed sets of rank r-1.
  I = [Set(0)] # I[r] is the family of independent sets of rank r-1.
  rank[0] = 0

  while E ∉ F[r]
    push!(F, Set())
    push!(I, Set())

    # Generate minimal closed sets for rank r+1.
    for y in F[r] # y is a closed set of rank r.
      t = E - y # The set of elements not in y.
      # Find all sets in F[r+1] that already contain y and remove excess elements from t.
      for x in F[r+1]
        if (x & y == y) t &= ~x end
      end
      # Insert y ∪ a for all a ∈ t.
      while t > 0
        x = y|(t&-t)
        insert_set!(x, F, r, rank)
        t &= ~x
      end
    end

    if E ∈ F[r+1]
      break
    end
    
    if r <= length(p)
      # Apply coarsening.
      pr = p[r]
      while pr > 0 && E ∉ F[r+1]
        A = rand(F[r+1])
        t = E-A
        one_element_added::Vector{T} = []
        while t > 0
          x = A|(t&-t)
          push!(one_element_added, x)
          t &= ~x
        end
        Acupa = rand(one_element_added)
        setdiff!(F[r+1], A)
        insert_set!(Acupa, F, r, rank)
        pr -= 1
      end
    end

    # Assign rank to sets and add independent ones to I.
    for m in F[r+1]
      mark!(m, I, r, rank)
    end

    # Next rank.
    r += 1
  end

  C = Set()
  k = 1
  while k <= E
    for i in 0:k-1 if rank[k+i] == rank[i]
      push!(C, T(k+i))
      unmark!(k+i, rank[i], rank, E)
    end end
    k += k
  end

  return KnuthMatroid{T}(n, F, I, C, rank)
end