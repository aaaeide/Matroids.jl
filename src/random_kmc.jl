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
  return random_bitwise_kmc(generate_covers_v2, superpose_v2!, n, p)
end

"""
Third implementation of random-KMC. This sorts the sets by size before superposing.
"""
function randomized_knuth_matroid_construction_v3(n, p)
  return random_bitwise_kmc(generate_covers_v2, sorted_superpose_v2!, n, p)
end