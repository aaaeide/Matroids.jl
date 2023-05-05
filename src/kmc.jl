include("types.jl")

"""
Knuth's matroid construction (1974). Generates a matroid in terms of its closed sets, given by the size of the universe n, a list of enlargements X and optionally the type T to use for set representation. T must have bitwidth >= n.
"""
function knuth_matroid(n, X, T=UInt16)
  @assert length(bitstring(T(n))) >= n 

  r = 1 # r is current rank +1 due to 1-indexing.
  F = [Set(0)]
  E = T(big"2"^n-1)
  rank = Dict{T, UInt8}(0 => 0)

  while E ∉ F[r]
    # Initialize F[r+1].
    push!(F, Set())

    # Setup add_set.
    add_callback = x -> rank[x] = r
    add_function = x -> add_set!(x, F, r, rank, add_callback)

    generate_covers!(F, r, E, add_function)
    
    # Perform enlargements.
    if r <= length(X) if X[r] !== nothing
      for x in X[r] add_function(x) end
    end end

    r += 1
  end

  return ClosedSetsMatroid{T}(n,r-1,F,rank, T)
end

"""
Randomized version of Knuth's matroid construction. Random matroids are generated via the method of random coarsening described in the 1974 paper. Accepts the size of the universe n, a list p = [p_1, p_2, ...], where p_i denotes the number of coarsenings to apply at rank i, and optionally the type T to use for set representation. T must have bitwidth >= n.
"""
function random_knuth_matroid(n, p, T=UInt16)::ClosedSetsMatroid{T}
  r = 1
  F::Vector{Set{T}} = [Set(T(0))]
  E::T = BigInt(2)^n-1
  rank = Dict{T, UInt8}(0=>0)

  while E ∉ F[r]
    # Initialize F[r+1].
    push!(F, Set())

    # Setup add_set.
    add_callback = x -> rank[x] = r
    add_function = x -> add_set!(x, F, r, rank, add_callback)

    generate_covers!(F, r, E, add_function)

    # Perform coarsening.
    if r <= length(p) coarsen!(F, r, E, p[r], add_function) end

    r += 1
  end

  return ClosedSetsMatroid{T}(n, r-1, F, rank, T) # , T)
end

"""
    Randomized version of Knuth's matroid construction, that also keeps track of independent sets. This entails keeping storing the whole power set of E in the rank table, and quickly becomes infeasible for values of n much larger than 16.
"""
function random_erect(n, p, T::Integer=UInt16)::FullMatroid{T}
  r = 1
  F::Vector{Set{T}} = [Set(T(0))]
  E::T = big"2"^n-1
  rank = Dict{T, UInt8}()

  # Populate rank table with 100+cardinality for all subsets of E.
  k=1; rank[0]=100;
  while (k<=E)
    for i in 0:k-1 rank[k+i] = rank[i]+1 end
    k=k+k;
  end

  F = [Set(0)]
  I = [Set(0)]
  rank[0] = 0

  while E ∉ F[r]
    # Initialize F[r+1] and I[r+1].
    push!(F, Set())
    push!(I, Set())

    # Setup add_set.
    add_callback = x -> mark_independent_subsets!(x, I, r, Base.count_ones(x), rank)
    add_function = x -> add_set!(x, F, r, rank, add_callback)
    
    generate_covers!(F, r, E, add_function)

    # Perform coarsening.
    if r <= length(p) coarsen!(F, r, E, p[r], add_function) end

    r += 1
  end

  return FullMatroid{T}(n, r-1, F, I, Set(), rank, T)
end

"""
Generates minimal closed sets for rank r+1 and inserts them into F[r+1], using the supplied insert_fn. This function should take one argument, the newly added set.
"""
function generate_covers!(F, r, E, insert_fn)
  for y in F[r]
    t = E - y
    # Find all sets in F[r+1] that already contain y and remove excess elements from t.
    for x in F[r+1]
      if (x & y == y) t &= ~x end
      if t == 0 break end
    end
    # Insert y ∪ a for all a ∈ t.
    while t > 0
      x = y|(t&-t)
      insert_fn(x)
      t &= ~x
    end
  end
end


function coarsen!(F, r, E, count, add_function)
  for _ in 1:count
    if E ∈ F[r+1] return end
    A = rand(F[r+1])
    t = E-A
    a = rand([1<<i for i in 0:length(bitstring(t)) if 1<<i & t != 0])
    setdiff!(F[r+1], A)
    add_function(A|a)
  end
end

function add_set!(x, F, r, rank, callback)
  for y in F[r+1]
    if haskey(rank, x&y) && rank[x&y]<r continue end

    # x ∩ y has rank > r, replace with x ∪ y.
    setdiff!(F[r+1], y)
    add_set!(x|y, F, r, rank, callback)
    return
  end

  push!(F[r+1], x)
  callback(x)
end

"""
Given a closed set x,
1. simply return if rank[x] < r (we've seen this already)
2. add it to I if |x| = r
3. recursively call this func on all x' ⊂ x st |x'| = |x| - 1
"""
function mark_independent_subsets!(x, I, r, c, rank)
  if haskey(rank, x) && rank[x] <= r return end
  if c == r push!(I[r+1], x) end
  rank[x] = r
  t = x
  while t != 0
    v = t&(t-1)
    mark_independent_subsets!(x-t+v, I, r, c-1, rank)
    t = v
  end
end
