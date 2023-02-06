# KnuthMatroid = Tuple(UInt16, Vector{Set{UInt16}})

function diff(A, B)
  A & ~B
end

function subseteq(A, B)
  A & B == A
end

function set_to_bits(set)
  sum(2^x for x in set)
end

function bits_to_set(bits)::Set{UInt16}
  return Set([i-1 for i in 1:16 if reverse(bitstring(bits))[i] == '1'])
end

function knuth_matroid_construction_v2(n, enlargements)#::KnuthMatroid
  # Step 1: Initialize.
  r = 1 # Julia is 1-indexed, subtract 1 to get current rank
  F::Vector{Set{UInt16}} = [Set(0)] # Start with the empty set.
  while true
    # Step 2: Generate covers.
    push!(F, generate_covers(F[r], n))

    # Step 3: Enlarge.
    if r <= length(enlargements) && enlargements[r] !== nothing
      F[r+1] = F[r+1] ∪ enlargements[r]
    end

    # Step 4: Superpose.

  end
end

"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the size of the universe.
"""
function generate_covers_v2(F_r, n)
  Set([A | 1 << i for A ∈ F_r for i in 0:n-1 if A & 1 << i === 0])
end

function superpose_v2!(F, F_prev)
  for A ∈ F
    println("OUTER\t", bits_to_set(A))
    for B ∈ F
      println("MIDDLE\t", bits_to_set(B))

      should_merge = true
      for C ∈ F_prev
        println("INNER\t", bits_to_set(C))
        if subseteq(A & B, C)
          println("\nMATCH FOUND! NOT MERGING\n")
          should_merge = false
        end

        readline()
      end

      if should_merge
        setdiff!(F, [A, B])
        push!(F, A | B)
      end
    end
  end

  return F
end