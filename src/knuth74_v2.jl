# KnuthMatroid = Tuple(UInt16, Vector{Set{UInt16}})

function diff(A, B)
  A & ~B
end

function subseteq(A, B)
  A & B == A
end

function random_element(A)
  rand([i-1 for i in 1:length(bitstring(A)) if reverse(bitstring(A))[i] === '1'])
end

function set_to_bits(set)
  UInt16(sum(2^x for x in set))
end

function bits_to_set(bits)
  return Set([i-1 for i in 1:16 if reverse(bitstring(bits))[i] == '1'])
end

function knuth_matroid_construction_v2(n, enlargements)
  # Step 1: Initialize.
  r = 1 # Julia is 1-indexed, subtract 1 to get current rank
  F::Vector{Set{UInt16}} = [Set(0)] # Start with the empty set.
  while true
    # Step 2: Generate covers.
    push!(F, generate_covers_v2(F[r], n))

    # Step 3: Enlarge.
    if r <= length(enlargements) && enlargements[r] !== nothing
      F[r+1] = F[r+1] ∪ enlargements[r]
    end

    # Step 4: Superpose.
    superpose_v2!(F[r+1], F[r])

    # Step 5: Test for completion.
    if 2^n-1 ∈ F[r+1]
      break
    end

    r += 1
  end

  return (n, F)
end

"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the size of the universe.
"""
function generate_covers_v2(F_r, n)
  Set([A | 1 << i for A ∈ F_r for i in 0:n-1 if A & 1 << i === 0])
end

function should_merge(A, B, F_prev)
  for C in F_prev
    # println("\n=====COMPARING=====")
    # println("OUTER\t", bits_to_set(A), "\t", bitstring(UInt16(A)))
    # println("MIDDLE\t", bits_to_set(B), "\t", bitstring(UInt16(B)))
    # println("INNER\t", bits_to_set(C), "\t", bitstring(UInt16(C)))
    
    if subseteq(A & B, C)
      # println("\nA ∩ B ⊆ C FOUND! NOT MERGING:")
      # println("A:\t", bits_to_set(A), "\t", bitstring(UInt16(A)))
      # println("B:\t", bits_to_set(B), "\t", bitstring(UInt16(B)))
      # println("----------------------------------------")
      # println("C:\t", bits_to_set(C), "\t", bitstring(UInt16(C)))
      # println("========================================")
      # readline()
      return false
    end

  end
  return true
end

function superpose_v2!(F, F_prev)
  As = collect(F)
  while length(As) !== 0
    @label next_A
    A = popfirst!(As)

    for B in setdiff(F, A)
      if should_merge(A, B, F_prev)
        insert!(As, 1, A | B)
        setdiff!(F, [A, B])
        push!(F, A | B)

        # println("\nA ∩ B ⊆ C NOT FOUND. MERGING ", bits_to_set(A), " AND ", bits_to_set(B))
        # println("REMOVING A ", bits_to_set(A), " AND B ", bits_to_set(B))
        # println("ADDING UNION ", bits_to_set(A | B))        
        # println("A in F == ", A in F)
        # println("B in F == ", B in F)        
        # readline()
        
        @goto next_A
      end
    end
  end

  return F
end


# function superpose_v2!(F, F_prev)
#   F_stk = sort!(collect(F), by = s -> length(bits_to_set(s)), rev = true)

#   while length(F_stk) != 0
#     A = popfirst!(F_stk)
#     # println("POPPING A: ", bits_to_set(A))

#     i = 1
#     while i <= length(F_stk)
#       B = F_stk[i]
#       if B === A
#         i += 1
#         continue
#       end

#       if should_merge(A, B, F_prev)
#         # Update stack.
#         deleteat!(F_stk, i) # Removes B from the stack. A is already popped.
#         if last(F_stk) != A | B
#           push!(F_stk, A | B)
#         end
        
#         # Update set.
#         setdiff!(F, [A, B])
#         push!(F, A | B)


        
#         # println("\nA ∩ B ⊆ C NOT FOUND. MERGING ", bits_to_set(A), " AND ", bits_to_set(B))
#         # println("REMOVING A ", bits_to_set(A), " AND B ", bits_to_set(B))
#         println("ADDING UNION ", bits_to_set(A | B))        
#         # println("A in F == ", A in F_stk)
#         # println("B in F == ", B in F_stk)        
#         # readline()
#       end
      
#       i+=1
#     end
#   end

#   return F
# # end

# function superpose_v2!(F, F_prev)
#   for A in F
    
#     for B in F
#       should_merge = true
#       for C in F_prev
#         if A & B & C == A & B # A ∩ B ⊆ C
#           should_merge = false
#         end
#       end

#       if should_merge
#         setdiff!(F, [A, B])
#         push!(F, A | B)
#       end
#     end
#   end
  
#   return F
# end


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
