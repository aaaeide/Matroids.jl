using ..Matroids
using  Allocations

## CHECKS ###########################################################

# General check for threshold-bases fairness measures (PROP, PROP1, PROPX, MMS),
# where each agent's bundle value should exceed some threshold. This threshold
# is computed with the supplied function thresh.
function check_thresh(V, A, thresh)
    for i in agents(V)
        if value(V, i, A) < thresh(V, i, A)
            return false
        end
    end

    return true
end

## ENVY-FREENESS ####################################################

# Matroid rank functions are binary, so if the bundle value 
# is > 0, we know we can remove a 1-valued item from the bundle.
value_1(V::MatroidRank, i, A) = max(value(V, i, A)-1, 0)

# EFX and EF1 is equivalent (the least positively valued item has the same 
# value as the highest valued item).
value_x(V::MatroidRank, i, A) = value_1(V, i, A)

# If S is independent for M_i, the least valued good in S has value 1, else 0.
value_x0(V::MatroidRank, i, A) =
    is_indep(V.matroids[i], A) ? value_1(V, i, A) : value(V, i, A)


## THRESHOLDS #######################################################

"""
    prop(V::MatroidRank, i, S)

The proportional share of agent `i`, that is, the value `i` places on the set 
of all goods, divided by the numder of agents.
"""
prop(V::MatroidRank, i, _) = rank(V.matroids[i])/na(V)

"""
    prop_1(V::Profile, i, S)

The proportional share of agent `i`, minus a highest-valued item not allocated
to `i`.
"""
prop_1(V::MatroidRank, i, A) = prop(V, i, A) - 1

"""
    prop_x(V::Profile, i, S)

The proportional share of agent `i`, minus a least positively-valued item not allocated to `i`.
"""
prop_x(V::MatroidRank, i, A) = prop_1(V, i, A)

"""
    prop_x0(V::Profile, i, S)

The proportional share of agent `i`, minus a least-valued (possibly 0-valued)
item not allocated to `i`.
"""
prop_x0(V::MatroidRank, i, A) = 
    is_closed(V.matroids[i], A) ? prop_1(V, i, A) : prop(V, i, A)

"""
    function mms_i(Mi, n)

Finds the maximin share of agent i for a fair allocation instance with matroid-rank valuations, given matroid Mi and the number of agents n.
"""
function mms_i(M_i, n)
  # An initial partition into independent subsets (subjectively so for i).
  A = matroid_partition_knuth_73([M_i for _ in 1:n])

  # Setup matrix D st D[j,k] v_i(A_j) - v_i(A_k) ∀ j,k ∈ [n].
  D = zeros(Int8, n, n)
  for j in 1:n for k in 1:n
    # v_i(A_p) = |A_p| since all sets in A are independent wrt M_i.
    D[j,k] = length(A[j]) - length(A[k])
  end end

  jk = argmax(D)
  while D[jk] > 1
    j,k = Tuple(jk)

    # By the augmentation property, ∃g ∈ A_j st A_k + g ∈ I_i.
    g = nothing
    for g´ ∈ setdiff(A[j], A[k]) if is_indep(M_i, A[k] ∪ g´)
      g = g´; break
    end end

    # Update A.
    setdiff!(A[j], g); union!(A[k], g)

    # Update D.
    for l in 1:n
      D[j, l] -= 1; D[l, j] += 1 # A_j is one smaller.
      D[k, l] += 1; D[l, k] -= 1 # A_k is one larger.
    end

    jk = argmax(D)
  end
  
  return minimum(length, A)
end