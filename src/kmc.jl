struct KnuthMatroid{T}
  n::Integer
  r::Integer
  F::Vector{Set{T}} # Closed sets by rank
  I::Vector{Set{T}} # Independent sets by rank
  C::Set{T} # Circuits
  rank::Dict{T, UInt8}
  Type::DataType
end

"""
Knuth's matroid construction (1974). Generates a matroid in terms of its closed sets, given by the size of the universe n, a list of enlargements X and optionally the type T to use for set representation. T must have bitwidth >= n.
"""
function knuth_matroid(n, X, T=UInt16)
  @assert length(bitstring(T(n))) >= n 

  r = 1
  F = [Set(0)]
  E = T(big"2"^n-1)
  rank = Dict{T, UInt8}(0 => 0)

  while E ∉ F[r]
    push!(F, Set())
    generate_covers_and_insert!(F, r, E, rank)
    if r <= length(X) enlarge!(F, r, X[r], rank) end
    r += 1
  end

  return KnuthMatroid{T}(n,r-1,F,[],Set(),rank, T)
end

"""
Randomized version of Knuth's matroid construction. Random matroids are generated via the method of random coarsening described in the 1974 paper. Accepts the size of the universe n, a list p = [p_1, p_2, ...], where p_i denotes the number of coarsenings to apply at rank i, and optionally the type T to use for set representation. T must have bitwidth >= n.
"""
function random_knuth_matroid(n, p, T=UInt16)::KnuthMatroid{T}
  r = 1
  F::Vector{Set{T}} = [Set(T(0))]
  E::T = BigInt(2)^n-1
  rank = Dict{T, UInt8}(0=>0)

  while E ∉ F[r]
    # Create empty set.
    push!(F, Set())
    generate_covers_and_insert!(F, r, E, rank)
    if r <= length(p) coarsen!(F, r, E, p[r], rank) end
    r += 1
  end

  return KnuthMatroid{T}(n, r-1, F, [], Set(), rank, T)
end

"""
Generates minimal closed sets for rank r+1 and inserts them into F[r+1].
"""
function generate_covers_and_insert!(F, r, E, rank)
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
      add_set!(x, F, r, rank)
      t &= ~x
    end
  end
end

function enlarge!(F, r, sets, rank)
  if sets === nothing return end
  for set in sets add_set!(set, F, r, rank) end
end

function coarsen!(F, r, E, count, rank)
  for _ in 1:count
    if E ∈ F[r+1] return end
    A = rand(F[r+1])
    t = E-A
    a = rand([1<<i for i in 0:length(bitstring(t)) if 1<<i & t != 0])
    setdiff!(F[r+1], A)
    add_set!(A|a, F, r, rank)
  end
end

function add_set!(x, F, r, rank)
  if x in F[r+1] return end
  for y in F[r+1]
    if haskey(rank, x&y) && rank[x&y]<r
      continue
    end

    # x ∩ y has rank > r, replace with x ∪ y.
    setdiff!(F[r+1], y)
    return add_set!(x|y, F, r, rank)
  end

  push!(F[r+1], x)
  rank[x] = r
end