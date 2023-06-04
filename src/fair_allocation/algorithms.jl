using Graphs
using Allocations

"""
    yankee_swap(V::MatroidRank)

Viswanathan and Zick's Yankee Swap algorithm (2022) for matroid-rank valuations. The agents are prioritized according to index.
"""
function yankee_swap(V::MatroidRank)
  A = falses(na(V)+1, ni(V))
  A[na(V)+1, :] .= 1 # The bundle of unallocated items.
  flag = falses(na(V))

  # Agent 0 has a corresponding zero matroid.
  D = exchange_graph([V.matroids..., ZeroMatroid(ni(V))], A)

  while false in flag
    # The agents whose bundle can still improve.
    T = [i for i in 1:na(V) if flag[i] == false]
    
    # Find the agents in T with minimim value.
    T_vals = [(i, value(V, i, A[i, :])) for i in T]
    min_val = minimum(last, T_vals)
    T´ = [i for (i, v) in T_vals if v == min_val]

    # The highest priority agent with minimum value.
    i = T´[1] 

    # The goods for which i has positive marginal value.
    F_i = [j for j in 1:ni(V) if Δ(V, i, A, j) == 1]

    # Find a shortest path from F_i to an unallocated good.
    A_0 = [j for j in 1:ni(V) if A[na(V)+1, j] == 1]
    transfer_path = find_shortest_path(D, F_i, A_0)

    # Transfer if path exists.
    if transfer_path !== nothing
      transfer!(A, i, transfer_path)
    else
      flag[i] = true
    end
  end

  alloc = Allocation(na(V), ni(V))
  for i in 1:na(V), j in 1:ni(V)
    if A[i,j] give!(alloc, i, j) end
  end
  
  return alloc
end



"""
    alloc_bciz21(V::MatroidRank)

The Envy-Induced Transfer algorithm, Algorithm 1 in Benabbou, Chakraborty, Igarashi and Zick (2021) computes a MAX-USW, EF1 
allocation.
"""
function alloc_bciz21(V::MatroidRank)
  n = na(V); m = ni(V)

  # Compute a clean, MAX-USW allocation.
  partition = matroid_partition_knuth73(V.matroids)
  A = Allocation(na(V), ni(V))
  for (i, bundle) in enumerate(partition)
    give!(A, i, bundle)
  end

  # D[i,j] holds v_i(A_i) - v_i(A_j).
  D = zeros(Int, n, n)
  for i in 1:n, j in 1:n
    D[i,j] = value(V, i, A; indep=true) - value(V, i, bundle(A, j))
  end

  # While there are agents i, j st i envies j more than 1...
  i,j = argmax(D) |> Tuple
  while D[i,j] > 1
    # Find item in A_j with marginal gain for i.
    o = findfirst(g -> Δ(V, A, i, g) == 1, collect(bundle(A, j)))

    # TROUBLE: o is nothing?? How can D[i,j] > 1 then?
    
    # Envy-induced transfer:
    deny!(A, j, o)
    give!(A, i, o)

    # Update D.
    for k in 1:n
      D[i, k] = value(V, i, A; indep=true) - value(V, i, bundle(A, k))
      D[k, i] = value(V, k, A; indep=true) - value(V, k, bundle(A, i))
      D[j, k] = value(V, j, A; indep=true) - value(V, j, bundle(A, k))
      D[k, j] = value(V, k, A; indep=true) - value(V, k, bundle(A, j))
    end

    i,j = argmax(D) |> Tuple
  end

  return A
end