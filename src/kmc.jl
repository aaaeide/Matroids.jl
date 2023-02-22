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

###################
# GENERATE COVERS #
###################

"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the ground set of elements E. Set-based.
"""
function generate_covers_v1(Fr, E)
  Set([A ∪ a for A ∈ Fr for a ∈ setdiff(E, A)])
end

"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the size of the universe. Bit-based.
"""
function generate_covers_v2(F_r, n)
  Set([A | 1 << i for A ∈ F_r for i in 0:n-1 if A & 1 << i === 0])
end

#############
# SUPERPOSE #
#############

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
  As = collect(F)
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

######################################
# KNUTH'S MATROID CONSTRUCTION (KMC) #
######################################
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

  while true
    # Step 2: Generate covers.
    push!(F, generate_covers_v1(F[r], E))

    # Step 3: Enlarge.
    if r <= length(enlargements) && enlargements[r] !== nothing
      F[r+1] = F[r+1] ∪ enlargements[r]
    end

    # Step 4: Superpose.
    superpose_v1!(F[r+1], F[r])

    # Step 5: Test for completion.
    if E ∈ F[r+1]
      break
    end

    r += 1
  end

  return (E, F)
end

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

function knuth_matroid_construction_v2(n, enlargements)
  return bitwise_kmc(generate_covers_v2, bitwise_superpose!, n, enlargements)
end

function knuth_matroid_construction_v3(n, enlargements)
  return bitwise_kmc(generate_covers_v2, sorted_bitwise_superpose!, n, enlargements)
end

"""
This is an attempt at a smarter implementation than directly following the setup from Knuth's 1974 article, instead inspired by Knuth's ERECTION.W implementation, wherein the superpose step is replaced by an insert operation that inserts new closed sets into the family of current rank one at a time, superposing on the fly.
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

struct KnuthMatroid{T}
  n::Integer
  F::Vector{Set{T}}
  rank::Dict{T, UInt8}
end

"""
In this implementation of KMC we are implementing a trick from Knuth's 
ERECTION.W, in which the covers are generated iteratively, a a rank table is used to keep track of set ranks, and sets are added one at a time.

The type of integer to use for representing sets as numbers can be supplied with the T argument. It defaults to UInt16, but needs to have at least as many bits as the size of the universe n.
"""
function knuth_matroid_construction_v6(n, enlargements, T=UInt16)::KnuthMatroid
  r = 1
  F = [Set(0)]
  mask = 2^n-1
  rank = Dict{T, UInt8}(0=>0)

  while mask ∉ F[r]
    # Create empty list.
    push!(F, Set())

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
        add_set!(x, F, r, rank)
        t &= ~x
      end
    end

    if r <= length(enlargements) && enlargements[r] !== nothing
      for set in enlargements[r]
        add_set!(set, F, r, rank)
      end
    end

    r += 1
  end

  return KnuthMatroid{T}(n,F,rank)
end

function add_set!(x, F, r, rank)
  for y in F[r+1]
    if haskey(rank, x&y) && rank[x&y]<r
      continue
    end

    if (x & y == x) return end # Not necessary?

    # x ∩ y has rank > r, replacex and y with x ∪ y.
    setdiff!(F[r+1], y)
    add_set!(x|y, F, r, rank)
    return
  end

  push!(F[r+1], x)
  rank[x] = r
end