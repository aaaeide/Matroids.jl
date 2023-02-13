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

      pr -= 1
    end

    # Step 5: Test for completion.
    if E ∈ F[r+1]
      break
    end

    r += 1
  end

  return (E, F)
end

"""
Second implementation of random-KMC.

n is the size of the universe.
p is a list (p_1, p_2, ...), where p_r is the number of coarsening steps to apply at rank r in the construction. The first entry of p should usually be 0, since adding closed sets of size > 1 at rank 1 is equivalent to shrinking E.

This uses the bit-based KMC methods.
"""
function randomized_knuth_matroid_construction_v2(n, p)
  # Step 1: Initialize.
  r = 1
  F = [Set(0x0000)]
  pr = 0
  E = UInt16(2^n-1)

  while true
    # Step 2: Generate covers.
    push!(F, generate_covers_v2(F[r], n))

    # Step 4: Superpose
    superpose_v2!(F[r+1], F[r])

    # Step 5: Test for completion.
    if E ∈ F[r+1]
      return (n, F)
    end

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
      superpose_v2!(F[r+1], F[r])
      
      # Step 5: Test for completion.
      if E ∈ F[r+1]
        return (n, F)
      end
      
      pr -= 1
    end


    r += 1
  end
end
