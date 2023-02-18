include("kmc.jl")

"""
First implementation of Knuth's random matroid construction through random "coarsening".

n is the size of the universe.
p is a list (p_1, p_2, ...), where p_r is the number of coarsening steps to apply at rank r in the construction. The first entry of p should usually be 0, since adding closed sets of size > 1 at rank 1 is equivalent to shrinking E.

This uses the Set-based KMC methods.
"""
function randomized_knuth_matroid_construction_v1(n, p)
  E = Set([i for i in range(0,n-1)])
  
  # Step 1: Initialize.
  r = 1
  F = [family([])]
  pr = 0
  
  while true
    # Step 2: Generate covers.
    push!(F, generate_covers_v1(F[r], E))
    
    # Step 4: Superpose.
    superpose_v1!(F[r+1], F[r])
    
    # Step 5: Test for completion.
    if E ∈ F[r+1]
      return (E, F)
    end
    
    if r <= length(p)
      pr = p[r]
    end
    
    while pr > 0
      # Random closed set in F_{r+1} and element in E ∖ A.
      A = rand(F[r+1])
      a = rand(setdiff(E, A))
      
      # Replace A with A ∪ {a}.
      F[r+1] = setdiff(F[r+1], A) ∪ Set([A ∪ a])
      
      # Superpose again to account for coarsening step.
      superpose_v1!(F[r+1], F[r])
      
      # Step 5: Test for completion.
      if E ∈ F[r+1]
        return (E, F)
      end
      
      pr -= 1
    end
    
    
    r += 1
  end
end


"""
Bitwise implementation of Knuth's approach to random matroid generation through a number of random "coarsening" steps. Supply the generate_covers and superpose methods to study the effects of different implementations of these.

n is the size of the universe.
p is a list (p_1, p_2, ...), where p_r is the number of coarsening steps to apply at rank r in the construction. The first entry of p should usually be 0, since adding closed sets of size > 1 at rank 1 is equivalent to shrinking E.
"""
function random_bitwise_kmc(generate_covers, superpose, n, p)
  # Initialize.
  r = 1
  pr = 0
  F = [Set(0)]
  E = 2^n - 1 # The set of all elements in E.

  while true
    # Generate covers.
    push!(F, generate_covers(F[r], n))

    # Superpose.
    superpose(F[r+1], F[r])

    # Test for completion.
    if E ∈ F[r+1]
      return (n, F)
    end

    # Apply coarsening.
    if r <= length(p)
      pr = p[r]
    end

    while pr > 0
      # Get random closed set A in F_{r+1} and element a in E ∖ A.
      A = rand(F[r+1])
      a = random_element(diff(E, A))

      # Replace A with A ∪ {a}.
      F[r+1] = setdiff(F[r+1], A) ∪ Set([A | a])

      # Superpose again to account for coarsening step.
      superpose(F[r+1], F[r])
      
      # Step 5: Test for completion.
      if E ∈ F[r+1]
        return (n, F)
      end
      
      pr -= 1
    end

    r += 1
  end
end

"""
Second implementation of random-KMC. This uses the bit-based KMC methods.
"""
function randomized_knuth_matroid_construction_v2(n, p)
  return random_bitwise_kmc(generate_covers_v2, bitwise_superpose!, n, p)
end

"""
Third implementation of random-KMC. This sorts the sets by size before superposing.
"""
function randomized_knuth_matroid_construction_v3(n, p)
  return random_bitwise_kmc(generate_covers_v2, sorted_bitwise_superpose!, n, p)
end

"""
Fourth implementation of random-KMC. This uses an iterative approach to superposition.
"""
function randomized_knuth_matroid_construction_v4(n, p)
  r = 1
  pr = 0
  F = [Set(0)]
  E = 2^n - 1 # The set of all elements in E.

  while true
    to_insert = collect(generate_covers_v2(F[r], n))

    # Apply coarsening to covers.
    if r <= length(p) && E ∉ to_insert # No need to coarsen if E is added.
      pr = p[r]
      while pr > 0
        A = rand(to_insert)
        a = random_element(E - A)
        to_insert = setdiff(to_insert, A) ∪ [A | a]
        pr -= 1
      end
    end

    # Superpose.
    push!(F, Set()) # Add F[r+1].
    while length(to_insert) > 0
      A = popfirst!(to_insert)
      push!(F[r+1], A)

      for B in setdiff(F[r+1], A)
        if should_merge(A, B, F[r])
          insert!(to_insert, 1, A | B)
          setdiff!(F[r+1], [A, B])
          push!(F[r+1], A | B)
        end
      end
    end

    if E ∈ F[r+1]
      return (n, F)
    end

    r += 1
  end
end

"""
Fifth implementation of random-KMC. This one uses a dictionary to keep track of previously seen sets.
"""
function randomized_knuth_matroid_construction_v5(n, p, T=UInt16)
  r = 1
  pr = 0
  F = [Set(0)]
  E = 2^n-1
  rank = Dict{T, UInt8}(0=>0) # The rank table maps from the representation of a set to its assigned rank.

  while true
    to_insert = collect(generate_covers_v2(F[r], n))

    # Apply coarsening to covers.
    if r <= length(p) && E ∉ to_insert # No need to coarsen if E is added.
      pr = p[r]
      while pr > 0
        A = rand(to_insert)
        a = random_element(E - A)
        to_insert = setdiff(to_insert, A) ∪ [A | a]
        pr -= 1
      end
    end

    # Superpose.
    push!(F, Set()) # Add F[r+1].
    while length(to_insert) > 0
      A = popfirst!(to_insert)
      push!(F[r+1], A)
      rank[A] = r

      for B in setdiff(F[r+1], A)
        if !haskey(rank, A&B) || rank[A&B] >= r
          # Update insert queue.
          insert!(to_insert, 1, A | B)

          # Update F[r+1].
          setdiff!(F[r+1], [A, B])
          push!(F[r+1], A | B)

          # Update rank table.
          rank[A|B] = r
          break
        end
      end
    end

    if E ∈ F[r+1]
      return (n, F)
    end

    r += 1
  end
end

"""
Sixth implementation of random-KMC, in which the covers are generated iteratively, a a rank table is used to keep track of set ranks, and sets are added one at a time.
"""

function randomized_knuth_matroid_construction_v6(n, p, T=UInt16)
  r = 1
  pr = 0
  F = [Set(0)]
  E = 2^n-1
  rank = Dict{T, UInt8}(0=>0)

  while E ∉ F[r]
    # Create empty set.
    push!(F, Set())
    
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
        add_set!(x, F, r, rank)
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
        a = random_element(E-A)
        setdiff!(F[r+1], A)
        add_set!(A|a, F, r, rank)
        pr -= 1
      end
    end

    r += 1
  end

  return (n, F)
end

function add_set!(x, F, r, rank)
  if x in F[r+1] return end
  for y in F[r+1]
    if haskey(rank, x&y) && rank[x&y]<r
      continue
    end

    # x ∩ y has rank > r, replace with x ∪ y.
    setdiff!(F[r+1], y)
    return add_set!(x|y, F, r, rank)
  end

  push!(F[r+1], x)
  rank[x] = r
end