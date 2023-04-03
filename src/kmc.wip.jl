##################
# VERSION ONE    #
##################

"""
  First implementation of Knuth's matroid construction algorithm (1974).
"""
function knuth_matroid_construction_v1(
  E,
  enlargements
)
  # Step 1: Initialize.
  r = 1
  F = [family([])]

  while E ∉ F[r]
    # Step 2: Generate covers.
    push!(F, generate_covers_v1(F[r], E))

    # Step 3: Enlarge.
    if r <= length(enlargements) && enlargements[r] !== nothing
      F[r+1] = F[r+1] ∪ enlargements[r]
    end

    # Step 4: Superpose.
    superpose_v1!(F[r+1], F[r])

    r += 1
  end

  return (E, F)
end

"""
First implementation of Knuth's random matroid construction through random "coarsening".

n is the size of the universe.
p is a list (p_1, p_2, ...), where p_r is the number of coarsening steps to apply at rank r in the construction. The first entry of p should usually be 0, since adding closed sets of size > 1 at rank 1 is equivalent to shrinking E.

This uses the Set-based KMC methods.
"""
function randomized_knuth_matroid_construction_v1(n, p, T)::KnuthMatroid{Set{Integer}}
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
      return KnuthMatroid{Set{Integer}}(n, F, [], Set(), Dict())
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
        return KnuthMatroid{Set{Integer}}(n, F, [], Set(), Dict())
      end
      
      pr -= 1
    end
    
    
    r += 1
  end
end

"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the ground set of elements E. Set-based.
"""
function generate_covers_v1(Fr, E)
  Set([A ∪ a for A ∈ Fr for a ∈ setdiff(E, A)])
end

"""
  If F_{r+1} contains any two sets A, B whose intersection A ∩ B is not contained in C for any C ∈ F_r, replace A, B ∈ F_{r+1} by the single set A ∪ B. Repeat this operation until A ∩ B ⊆ C for some C ∈ F_r whenever A and B are distinct members of F_{r+1}.

  F and F_old should be Family: A Set of Sets of some type.

  This is the first, Set-based implementation of this method.
"""
function superpose_v1!(F, F_old)
  for A ∈ F
    for B ∈ F
      should_merge = true
      for C ∈ F_old
        if A ∩ B ⊆ C
          should_merge = false
        end
      end

      if should_merge
        setdiff!(F, [A, B])
        push!(F, A ∪ B)
      end
    end
  end

  return F
end






##################
# VERSION TWO    #
##################

"""
Knuth's matroid construction (KMC) algorithm. In this version, the sets are represented using bits.
"""
function bitwise_kmc(generate_covers, superpose!, n, enlargements)
  # Step 1: Initialize.
  r = 1 # Julia is 1-indexed, subtract 1 to get current rank
  F = [Set(0)] # Start with the empty set.

  while 2^n-1 ∉ F[r]
    # Step 2: Generate covers.
    push!(F, generate_covers(F[r], n))

    # Step 3: Enlarge.
    if r <= length(enlargements) && enlargements[r] !== nothing
      F[r+1] = F[r+1] ∪ enlargements[r]
    end

    # Step 4: Superpose.
    superpose!(F[r+1], F[r])
    r += 1
  end

  return (n, F)
end


"""
Bitwise implementation of Knuth's approach to random matroid generation through a number of random "coarsening" steps. Supply the generate_covers and superpose methods to study the effects of different implementations of these.

n is the size of the universe.
p is a list (p_1, p_2, ...), where p_r is the number of coarsening steps to apply at rank r in the construction. The first entry of p should usually be 0, since adding closed sets of size > 1 at rank 1 is equivalent to shrinking E.
"""
function random_bitwise_kmc(generate_covers, superpose, n, p)::KnuthMatroid{Any}
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
      return KnuthMatroid{Any}(n, F, [], Set(), Dict())
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
        return KnuthMatroid{Any}(n, F, [], Set(), Dict())
      end
      
      pr -= 1
    end

    r += 1
  end
end

"""
Second implementation of random-KMC. This uses the bit-based KMC methods.
"""
function randomized_knuth_matroid_construction_v2(n, p, T=UInt16)
  return random_bitwise_kmc(generate_covers_v2, bitwise_superpose!, n, p)
end

function knuth_matroid_construction_v2(n, enlargements)
  return bitwise_kmc(generate_covers_v2, bitwise_superpose!, n, enlargements)
end

"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the size of the universe. Bit-based.
"""
function generate_covers_v2(F_r, n)
  Set([A | 1 << i for A ∈ F_r for i in 0:n-1 if A & 1 << i === 0])
end

"""
Returns whether the intersection of A and B is contained within
some C in F_prev.
"""
function should_merge(A, B, F_prev)
  for C in F_prev
    if subseteq(A & B, C)
      return false
    end
  end
  return true
end

"""
If F contains any two sets A, B whose intersection A ∩ B is not contained in C for any C ∈ F_prev, replace A, B ∈ F with the single set A ∪ B. Repeat this operation until A ∩ B ⊆ C for some C ∈ F_prev whenever A and B are distinct members of F.

This implementation represents the sets using bits.
"""
function bitwise_superpose!(F, F_prev)
  i = 0
  As = copy(F)
  while length(As) !== 0
    A = pop!(As)

    for B in setdiff(F, A)
      i += 1
      if should_merge(A, B, F_prev)
        push!(As, A | B)
        setdiff!(F, [A, B])
        push!(F, A | B)
        break
      end
    end
  end

  return F
end







##################
# VERSION THREE  #
##################

function knuth_matroid_construction_v3(n, enlargements)
  return bitwise_kmc(generate_covers_v2, sorted_bitwise_superpose!, n, enlargements)
end

"""
Third implementation of random-KMC. This sorts the sets by size before superposing.
"""
function randomized_knuth_matroid_construction_v3(n, p, T=UInt16)
  return random_bitwise_kmc(generate_covers_v2, sorted_bitwise_superpose!, n, p)
end

function sorted_bitwise_superpose!(F, F_prev)
  As = sort!(collect(F), by = s -> length(bits_to_set(s)))
  while length(As) !== 0
    A = popfirst!(As)

    for B in setdiff(F, A)
      if should_merge(A, B, F_prev)
        insert!(As, 1, A | B)
        setdiff!(F, [A, B])
        push!(F, A | B)
        break
      end
    end
  end

  return F
end






##################
# VERSION FOUR   #
##################

"""
This is an attempt at a smarter implementation than directly following the setup from Knuth's 1974 article. The superpose step is replaced by an insert operation that inserts new closed sets into the family of current rank one at a time, superposing on the fly.
"""
function knuth_matroid_construction_v4(n, enlargements)
  r = 1
  F = [Set(0)]
  
  while 2^n-1 ∉ F[r]
    to_insert = collect(generate_covers_v2(F[r], n))
    if r <= length(enlargements) && enlargements[r] !== nothing
      append!(to_insert, enlargements[r])
    end

    push!(F, Set())  # Add F[r+1]
    while length(to_insert) > 0
      A = popfirst!(to_insert)
      push!(F[r+1], A)

      for B in setdiff(F[r+1], A)
        if should_merge(A, B, F[r])
          insert!(to_insert, 1, A | B)
          setdiff!(F[r+1], [A, B])
          push!(F[r+1], A | B)
          break
        end
      end
    end

    r += 1
  end

  return (n, F)
end

"""
Fourth implementation of random-KMC. The superpose step is replaced by an insert operation that inserts new closed sets into the family of current rank one at a time, superposing on the fly.
"""
function randomized_knuth_matroid_construction_v4(n, p, T=UInt16)::KnuthMatroid{T}
  r = 1
  pr = 0
  F = [Set(0)]
  E = 2^n - 1 # The set of all elements in E.

  while true
    to_insert = generate_covers_v2(F[r], n)

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
      A = pop!(to_insert)
      push!(F[r+1], A)

      for B in setdiff(F[r+1], A)
        if should_merge(A, B, F[r])
          push!(to_insert, A | B)
          setdiff!(F[r+1], [A, B])
          push!(F[r+1], A | B)
        end
      end
    end

    if E ∈ F[r+1]
      return KnuthMatroid{T}(n, F, [], Set(), Dict())
    end

    r += 1
  end
end






##################
# VERSION FIVE   #
##################

"""
In a sense memoized with a rank table.
"""
function knuth_matroid_construction_v5(n, enlargements, T=UInt16)
  r = 1
  F = [Set(0)]
  rank = Dict{T, UInt8}(0=>0) # The rank table maps from the representation of a set to its assigned rank.
  
  while 2^n-1 ∉ F[r]
    to_insert = collect(generate_covers_v2(F[r], n))
    if r <= length(enlargements) && enlargements[r] !== nothing
      append!(to_insert, enlargements[r])
    end

    push!(F, Set())  # Add F[r+1]
    while length(to_insert) > 0
      A = popfirst!(to_insert)
      push!(F[r+1], A)
      rank[A] = r # Julia is 1-indexed, so F[r+1] is the family of rank r.

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

    r += 1
  end

  return (n, F)
end


"""
Fifth implementation of random-KMC. This one uses a dictionary to keep track of previously seen sets.
"""
function randomized_knuth_matroid_construction_v5(n, p, T=UInt16)::KnuthMatroid{T}
  r = 1
  pr = 0
  F::Vector{Set{T}} = [Set(T(0))]
  E = 2^n-1
  rank = Dict{T, UInt8}(0=>0) # The rank table maps from the representation of a set to its assigned rank.

  while true
    to_insert = generate_covers_v2(F[r], n)

    # Apply coarsening to covers.
    if r <= length(p)
      pr = p[r]
      while length(to_insert) > 0 && pr > 0 && E ∉ to_insert # No need to coarsen if E is added.
        A = rand(to_insert)
        a = random_element(E - A)
        to_insert = setdiff(to_insert, A) ∪ [A | a]
        pr -= 1
      end
    end

    # Superpose.
    push!(F, Set()) # Add F[r+1].
    while length(to_insert) > 0
      A = pop!(to_insert)
      push!(F[r+1], A)
      rank[A] = r

      for B in setdiff(F[r+1], A)
        if !haskey(rank, A&B) || rank[A&B] >= r
          # Update insert queue.
          push!(to_insert, A | B)

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
      return KnuthMatroid{T}(n, F, [], Set(), rank)
    end

    r += 1
  end
end






##################
# HELPER METHODS #
##################


"""
Generate a family - a set of sets, each inner set the result of calling
Set(x), where x is an element in the supplied Vector.
"""
function family(set_contents::Vector{Any})
  Set([Set([])])
end

function family(set_contents)
  Set([Set(e) for e in set_contents])
end

function pretty_print_family(F)
  for i in 1:(length(F)-1)
    println("\n\nF$i")
    print("{")
    for s in F[i]
      for e in s
        print(e)
      end
      print(", ")
    end
    print("}\n\n")
  end
end

function diff(A, B)
  A & ~B
end

function subseteq(A, B)
  A & B == A
end

function random_element(A)
  rand([i-1 for (i,c) in enumerate(reverse(bitstring(A))) if c == '1'])
end

function set_to_bits(set)
  sum(2^x for x in set)
end

function bits_to_set(bits)
  Set(i-1 for (i, c) in enumerate(reverse(bitstring(bits))) if c == '1')
end

