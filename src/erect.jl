"""
An improved version of KMC we are also finding the independent sets and circuits of the matroid during generation.

This version assigns the Hamming weight of all subsets of E upfront. This is infeasible for values of n much larger than 16.
"""
function erect_v1(n, enlargements)::KnuthMatroid{UInt16}
  # Initialize.
  r = 1
  mask = 2^n-1
  rank = Dict{UInt16, UInt8}()
  
  # Populate rank table with 100+cardinality for all subsets of E.
  k=1; rank[0]=100;
  while (k<=mask)
    for i in 0:k-1 rank[k+i] = rank[i]+1 end
    k=k+k;
  end
  
  F = [Set(0)] # F[r] is the family of closed sets of rank r-1.
  I = [Set(0)] # I[r] is the family of independent sets of rank r-1.
  rank[0] = 0

  while mask ∉ F[r]
    push!(F, Set())
    push!(I, Set())

    # Generate minimal closed sets for rank r+1.
    for y in F[r] # y is a closed set of rank r.
      t = mask - y # The set of elements not in y.
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

    # Enlarge (if any).
    if r <= length(enlargements) && enlargements[r] !== nothing
      for set in enlargements[r]
        insert_set!(set, F, r, rank)
      end
    end 
    
    # Assign rank to sets and add independent ones to I.
    for m in F[r+1]
      mark!(m, I, r, rank)
    end
    
    # Next rank.
    r += 1
  end

  C = Set() # C is the set of circuits (minimal dependent sets) for M.
  k = 1
  while k <= mask
    for i in 0:k-1 if rank[k+i] == rank[i]
      push!(C, UInt16(k+i))
      unmark!(k+i, rank[i]+101, rank, mask)
    end end
    k += k
  end

  return KnuthMatroid{UInt16}(n,F,I,C,rank)
end

"""
Inserts set x into F[r], but augments x if it is necessary to ensure no two sets in F[r] have an intersection of rank greater than r.
"""
function insert_set!(x, F, r, rank)
  for y in F[r+1] # +1 since Julia is 1-indexed.
    if rank[x&y] <= r
      continue
    end

    # x ∩ y has rank > r, replace x and y with x ∪ y.
    setdiff!(F[r+1], y)
    insert_set!(x|y, F, r, rank)
    return
  end

  push!(F[r+1], x)
end

"""
Given a closed set m, sets rank[m']=r for all subsets m' ⊆ m whose rank is not already ≤ r, and adds m' to I if it is independent (that is, if its rank equals its cardinality).
"""
function mark!(m, I, r, rank)
  if haskey(rank, m) && rank[m] <= r
    return
  end
  if rank[m] == 100+r push!(I[r+1], m) end
  rank[m] = r
  t = m
  while t != 0
    v = t&(t-1)
    mark!(m-t+v, I, r, rank)
    t = v
  end
end

function unmark!(m, card, rank, mask)
  if rank[m] < 100
    rank[m] = card
    t = mask-m
    while t != 0
      v = t&(t-1)
      unmark!(m+t-v, card+1, rank, mask)
      t=v
    end
  end
end

function set_to_string(t)
  str = ""
  for j in 0:16 if t&(1<<j)!=0 str *= string(j) end end
  return str * " "
end

"""
Generates minimal closed sets for rank r+1 and inserts to F[r+1] using supplied insert function.
  """
  function generate_covers_and_insert(F, r, E, rank, ins!)
    for y in F[r] # y is a closed set of rank r.
    t = E - y   # The set of elements not in y.
    # Find all sets in F[r+1] that already contain y and remove excess elements from t.
    for x in F[r+1]
      if (x&y == y) t &= ~x end # x subseteq y
    end
    # Callback y cup a for all a in t.
    while t > 0
      x = y|(t&-t)
      ins!(x, F, r, rank)
      t &= ~x
    end
  end
end

"""
Inserts set x into F[r+1], but augments x if it is necessary to ensure no two
sets in F[r] have an intersection of rank greater than r. This ensures that 
F[r] only consists of closed sets (maximal dependent sets).

Adds closed sets to rank table.
"""
function insert_set_v2!(x, F, r, rank)
  for y in F[r+1]
    # If x&y is in the rank table, it has rank <= r.
    if haskey(rank, x&y) && rank[x&y] < r continue end

    # x&y has rank > r (not seen yet), replace x and y with x|y.
    setdiff!(F[r+1], y)
    insert_set_v2!(x|y, F, r, rank)
    return
  end

  push!(F[r+1], x)
  rank[x] = r
end

function enlarge!(enlargements, F, r, rank, ins!)
  if r <= length(enlargements) && enlargements[r] !== nothing
    for set in enlargements[r] ins!(set, F, r, rank) end
  end
end

"""
Given a closed set m,
1. simply return if rank[m] < r (we've seen this already)
2. add it to I if |m| = r
3. recursively call this func on all m' ⊂ m st |m'| = |m| - 1
"""
function mark_independent_sets!(m, I, r, rank)
  if haskey(rank, m) && rank[m] < r return end
  if Base.count_ones(m) == r push!(I[r+1], m) end
  rank[m] = r

  t = m
  while t != 0
    v = t&(t-1)
    mark_independent_sets!(m-t+v, I, r, rank)
    t = v
  end
end


function erect_v2(n, enlargements, T=UInt16)::KnuthMatroid{T}
  r = 1
  E = big"2"^n-1
  rank = Dict{T, UInt8}() # Keeps track of closed sets' ranks.

  F = [Set(0)]
  I = [Set(0)]
  rank[0] = 0

  while E ∉ F[r]
    push!(F, Set())
    push!(I, Set())

    generate_covers_and_insert(F, r, E, rank, insert_set_v2!)
    enlarge!(enlargements, F, r, rank, insert_set_v2!)

    for m in F[r+1] mark_independent_sets!(m, I, r, rank) end

    r += 1
  end

  return KnuthMatroid{T}(n, F, I, Set(), rank)
end