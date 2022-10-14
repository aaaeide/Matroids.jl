using DataStructures

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
Generate the set F_{r+1} of all "covers" of the sets in F_r.
"""
function generate_covers(Fr::Family{<:Any}, E::Set{<:Any})::Family{<:Any}
  Set([A ∪ a for A ∈ Fr for a ∈ setdiff(E, A)])
end

# TODO: Make work with <:Any
function enlarge(family::Family{<:Any}, to_add::Vector{Set{Int64}})::Family{<:Any}
  family ∪ to_add
end

"""
  If F_{r+1} contains any two sets A, B whose intersection A n B is not contained in C for any C in F_r, replace A and B in F_{r+1} by the single set A U B. Repeat this operation until A n B subseteq C for some C in F_r whenever A and B are distinct members of F_{r+1}.
"""
function superpose(f_cur::Family{<:Any}, f_old::Family{<:Any})::Family{<:Any}
  f_new = Family{<:Any}()
  stack = Stack{Set{<:Any}}()
  for set ∈ f_cur
    push!(stack, set)
  end

  while !isempty(stack)
    A, B = pop!(stack), pop!(stack)
    for C ∈ f_old
      if A ∩ B ∈ C
        push!(f_new, A)
        push!(f_new, B)
        continue
      end
    end

    push!(stack, A ∪ B)
  end

  return f_new
end

export knuth_matroid_construction, generate_covers