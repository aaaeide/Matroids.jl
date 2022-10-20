Family{T} = Set{Set{T}}

"""
  First implementation of Knuth's matroid construction algorithm (1974).
"""
function knuth_matroid_construction(
  E::Set{Integer}, 
  enlargements:: Vector{Vector{Set{<:Any}}}
)
  # Step 1: Initialize.
  r = 0 
  F::Array{Family{Integer}} = [Set(Set())]

  # Step 2: Generate covers.
  
end

"""
Generate the set F_{r+1} of all "covers" of the sets in F_r, given the ground set of elements E.
"""
function generate_covers(Fr, E)::Family{<:Any}
  Set([A ∪ a for A ∈ Fr for a ∈ setdiff(E, A)])
end

"""
Enlarge a family by adding some additional sets of elements.
"""
function enlarge(family, to_add)::Family{<:Any}
  family ∪ to_add
end

"""
  If F_{r+1} contains any two sets A, B whose intersection A ∩ B is not contained in C for any C ∈ F_r, replace A, B ∈ F_{r+1} by the single set A ∪ B. Repeat this operation until A ∩ B ⊆ C for some C ∈ F_r whenever A and B are distinct members of F_{r+1}.

  f_cur and f_old should be Family: A Set of Sets of some type
"""
function superpose(f_cur, f_old)::Family{<:Any}
  f_new = Set(f_cur)

  for A ∈ f_new
    for B ∈ f_new
      should_merge = true
      for C ∈ f_old
        if A ∩ B ⊆ C
          should_merge = false
        end
      end

      if should_merge
        setdiff!(f_new, [A, B])
        push!(f_new, A ∪ B)
      end
    end
  end

  return f_new
end

