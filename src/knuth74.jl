Family{T} = Set{Set{T}}
KnuthMatroid{T} = Tuple{Set{T}, Vector{Family{T}}}


"""
Generate a family - a set of sets, each inner set the result of calling
Set(x), where x is an element in the supplied Vector.
"""
function family(set_contents::Vector{Any})::Family{Any}
  Set([Set([])])
end

function family(set_contents)::Family{<:Any}
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

"""
  First implementation of Knuth's matroid construction algorithm (1974).
"""
function knuth_matroid_construction(
  E,
  enlargements
)::KnuthMatroid{Any}
  # Step 1: Initialize.
  r = 1
  F = [family([])]

  while true
    # Step 2: Generate covers.
    push!(F, generate_covers(F[r], E))

    # Step 3: Enlarge.
    if r <= length(enlargements) && enlargements[r] !== nothing
      F[r+1] = F[r+1] ∪ enlargements[r]
    end

    # Step 4: Superpose.
    superpose!(F[r+1], F[r])

    # Step 5: Test for completion.
    if E ∈ F[r+1]
      break
    end

    r += 1
  end

  return (E, F)
end

"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the ground set of elements E.
"""
function generate_covers(Fr, E)::Family{<:Any}
  Set([A ∪ a for A ∈ Fr for a ∈ setdiff(E, A)])
end

"""
  If F_{r+1} contains any two sets A, B whose intersection A ∩ B is not contained in C for any C ∈ F_r, replace A, B ∈ F_{r+1} by the single set A ∪ B. Repeat this operation until A ∩ B ⊆ C for some C ∈ F_r whenever A and B are distinct members of F_{r+1}.

  F and F_old should be Family: A Set of Sets of some type
"""
function superpose!(F, F_old)::Family{<:Any}
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
Find the bases of a matroid given an array F of its closed sets grouped by rank,  F = [F_0, F_1, ..., F_n], where n is the rank of the matroid.
"""
function bases(F)#::Family{<:Any}
  rank = length(F) - 1 # Julia is 1-indexed, so F[1] = F_0
  elements = pop!(F[length(F)])
  min_closed_sets_of_2nd_highest_rank = 
    [s for s in F[length(F) - 1] if length(s) == rank - 1]

  bases = Set()
  for s in min_closed_sets_of_2nd_highest_rank
    for e in elements
      b = union(s, e)
      if length(b) == rank
        push!(bases, b)
      end
    end
  end

  for b in sort!([join([string(x) for x in sort(collect(b))]) for b in bases])
    println(b)
  end

  return bases
end