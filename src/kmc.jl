include("types.jl")

function VERBOSE(thresh) thresh >= 100 end
function bst(x) bitstring(x)[1:end] end

"""
Knuth's matroid construction (1974). Generates a matroid in terms of its closed sets, given by the size of the universe n, a list of enlargements X and optionally the type T to use for set representation. T must have bitwidth >= n.
"""
function knuth_matroid(n, X, T=UInt16)
  @assert length(bitstring(T(n))) >= n

  r = 1 # r is current rank +1 due to 1-indexing.
  F = [Set(0)]
  E = T(big"2"^n-1)
  rank = Dict{T, Integer}(0 => 0)

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
function random_knuth_matroid(n, p, T=UInt16, OVERRIDE=[])::ClosedSetsMatroid{T}
  r = 1
  F::Vector{Set{T}} = [Set(T(0))]
  E::T = BigInt(2)^n-1
  rank = Dict{T, Integer}(0=>0)

  while E ∉ F[r]
    # Initialize F[r+1].
    push!(F, Set())

    # Setup add_set.
    add_callback = x -> rank[x] = r
    add_function = x -> add_set!(x, F, r, rank, add_callback)

    generate_covers!(F, r, E, add_function)

    # Perform coarsening.
    if length(OVERRIDE) == 0
      if r <= length(p) coarsen!(F, r, E, p[r], add_function) end
    else
      if r <= length(OVERRIDE) coarsen_exact!(F, r, OVERRIDE[r], add_function) end
    end

    r += 1
    VERBOSE(1) && readline()
  end

  return ClosedSetsMatroid{T}(n, r-1, F, rank, T)
end

"""
    Randomized version of Knuth's matroid construction, that also keeps track of independent sets. This entails keeping storing the whole power set of E in the rank table, and quickly becomes infeasible for values of n much larger than 16.
"""
function random_erect(n, p, T=UInt16, OVERRIDE=[])::FullMatroid{T}
  r = 1
  E::T = big"2"^n-1
  rank = Dict{T, Integer}(0=>0)

  F::Vector{Set{T}} = [Set(T(0))]
  I::Vector{Set{T}} = [Set(T(0))]

  coarsenings = []

  while E ∉ F[r]
    # Initialize F[r+1] and I[r+1].
    push!(F, Set())
    push!(I, Set())

    # Setup add_set.
    add_callback = x -> mark_independent_subsets!(x, I, r, Base.count_ones(x), rank)
    add_function = x -> add_set!(x, F, r, rank, add_callback)
  
    generate_covers!(F, r, E, add_function)

    # Perform coarsening.
    if length(OVERRIDE) == 0
      if r <= length(p) coarsen!(F, r, E, p[r], add_function, coarsenings) end
    else
      if r <= length(OVERRIDE) coarsen_exact!(F, r, OVERRIDE[r], add_function) end
    end

    r += 1
    VERBOSE(1) && readline()
  end

  return FullMatroid{T}(n, r-1, F, I, Set(), rank, T)
end

function coarsen_exact!(F, r, Aas, add_fn)
  for (A, a) in Aas
    VERBOSE(2) && println("\n********  coarsening (exact)  ********")
    setdiff!(F[r+1], A)
    add_fn(A|a)
  end
end

"""
Generates minimal closed sets for rank r+1 and inserts them into F[r+1], using the supplied insert_fn. This function should take one argument, the newly added set.
"""
function generate_covers!(F, r, E, insert_fn)
  for y in F[r]
    VERBOSE(2) && println("\n=== generating covers for $(bst(y)) ===")
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


function coarsen!(F, r, E, count, add_function, log=[])
  entry = []
  for _ in 1:count
    if E ∈ F[r+1] return end
    A = rand(F[r+1])
    t = E-A
    a = rand([1<<i for i in 0:length(bst(t)) if 1<<i & t != 0])
    a = convert(typeof(A), a)
    setdiff!(F[r+1], A)

    VERBOSE(2) && println("\n********      coarsening      ********")
    add_function(A|a)

    push!(entry, (A, a))
  end
  push!(log, entry)
end

function add_set!(x, F, r, rank, callback)
  VERBOSE(1) && readline()
  VERBOSE(2) && print(" trying to add $(bst(x)) to rank $r")
  VERBOSE(2) && if Base.count_ones(x) > r print(" (REDUNDANCE)\n") else print("\n") end
  for y in F[r+1]
  
    """
    When we are comparing a set X which we are investigating whether is a closed set, with a set Y which we know (based on the sets that have been added thus far) is a closed set, we can encounter three scenarios:

    1. X ∩ Y is a closed set (it is in our rank table) of rank < r
      - Move on. If this is the case for all Y ∈ F[r+1] we add X to F[r+1].
    2. X ∩ Y is a closed set (it is in our rank table) of rank == r
      - X and Y are not closed sets. Remove Y from F[r+1] and call this
        function on X ∪ Y.
    3. We have not seen X ∩ Y before (this happens when we have closed sets of
       lower rank but similar cardinality).
      a. If |X ∩ Y| < r, we know that the rank of X ∩ Y is < r. Move on.
      b. If |X ∩ Y| >= r, we need to see if X ∩ Y ⊆ Z for some Z of lower rank.
         If not, remove Y from F[r+1] and call this function on X ∪ Y.

    """

    if haskey(rank, x&y) && rank[x&y]<r
      VERBOSE(0) && println("\ttable ok ($(rank[x&y])): $(bst(y))")
      continue
    end

    if !haskey(rank, x&y)
      if Base.count_ones(x&y) < r
        VERBOSE(0) && println("\tcardinality ok: $(bst(y))")
        continue
      else
        VERBOSE(4) && println("\tBAD CARDINALITY: $(bst(x&y)) - CHECKING...")
      
        r´ = check_rank(x&y, r, F)
        if r´ !== false
          rank[x&y] = r´
          continue
        end
      end
    end

    VERBOSE(2) && println("\t$(bst(x)) ∩ $(bst(y)) = $(bst(x&y)) has rank == $r")
    VERBOSE(2) && println("\treplacing with $(bst(x|y))")

    # x ∩ y has rank > r, replace with x ∪ y.
    setdiff!(F[r+1], y)
    add_set!(x|y, F, r, rank, callback)
    return
  end

  VERBOSE(3) && println("\tadding $(bst(x)) to rank $r")

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

function check_rank(v, r, F)
  for (i, Fi) in enumerate(F[1:r])
    for z ∈ Fi
      if v&z == v
        VERBOSE(4) && println("\tCHECK OKAY")
        return i-1
      end
    end
  end
  VERBOSE(4) && println("\tCHECK NOT OKAY. MERGE!")
  return false
end