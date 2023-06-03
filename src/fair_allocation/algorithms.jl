"""
    yankee_swap(V::MatroidRank)

Viswanathan and Zick's Yankee Swap algorithm (2022) for matroid-rank valuations. The agents are prioritized according to index.
"""
function yankee_swap(V::MatroidRank)
  A = falses(na(V)+1, ni(V))
  A[na(V)+1, :] .= 1 # The bundle of unallocated items.
  flag = falses(na(V))

  while false in flag
    T = [i for i in 1:na(V) if flag[i] == false]
    T_vals = [(i, value(V, i, A[i, :])) for i in T]
    min_val = minimum(last, T_vals)
    T´ = [i for (i, v) in T_vals if v == min_val]
    i = T´[1] # The highest priority agent.


  end
end