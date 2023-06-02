include("types.jl")
include("utils.jl")

ground_set(M::ClosedSetsMatroid) = bits_to_set(2^M.n-1)
ground_set(M::FullMatroid) = bits_to_set(2^M.n-1)
ground_set(M::UniformMatroid) = Set([1:M.n])



"""
    is_indep(M::Matroid, S::Set)

Independence oracle. Determines whether S is independent in M.
"""
function is_indep end

"""
    is_indep(M::ClosedSetsMatroid, S::Integer)

Determines whether a given set S is independent in the matroid M, given by the closed sets of M grouped by rank. Uses (I.1) in Greene (1989).
"""
function is_indep(M::ClosedSetsMatroid, S::Integer)
  t = Base.count_ones(S)

  if t > length(M.F) return false end

  for F in M.F[t]
    if Base.count_ones(S&F) > t-1 return false end
  end

  return true
end

is_indep(M::ClosedSetsMatroid, S) = is_indep(M, set_to_bits(S))


function is_indep(M::FullMatroid, S::Integer)
  card = Base.count_ones(S)
  if card + 1 > length(M.I) return false end
  return S in M.I[card+1]
end

is_indep(M::FullMatroid, S) = is_indep(M, set_to_bits(S))


function is_indep(M::UniformMatroid, S::Integer)
  return Base.count_ones(S) <= M.r
end

is_indep(M::UniformMatroid, S) = is_indep(M, set_to_bits(S))



"""
    is_circuit(M::Matroid, S)

Determines whether S is a circuit in M (ie. rank(S) = |S|-1).
"""
function is_circuit end


"""
    function is_circuit(M::ClosedSetsMatroid, S::Integer)

Determines whether a given set S is a circuit in the matroid M, given by the closed sets of M. Uses (C.1) and (C.2) in Greene (1989).
"""
function is_circuit(M::ClosedSetsMatroid, S::Integer)
  t = Base.count_ones(S)

  for F in M.F[t] # (C.1) S ⊆ F for some F ∈ F_{t-1}.
    if S&F==S @goto C2 end
  end
  return false

  @label C2
  for F in M.F[t-1] # (C.2) |S ∩ F| ≤ r(F) for all F ∈ F_{t-2}.
    if Base.count_ones(S&F) > t-2 return false end
  end

  return true
end


is_circuit(M::FullMatroid, S::Integer) = S in M.C
is_circuit(M::FullMatroid, S) = is_circuit(M, set_to_bits(S))


is_circuit(M::UniformMatroid, S::Integer) = r == Base.count_ones(S)-1 <= M.n
is_circuit(M::UniformMatroid, S) = is_circuit(M, set_to_bits(S))



"""
    minimal_spanning_subset(M::Matroid, S)

Finds a minimal spanning subset of S in M. If S is the ground set of M, this produces a basis of M.
"""
function minimal_spanning_subset end


"""
    function minimal_spanning_subsets(M::ClosedSetsMatroid, A::Integer)

Algorithm 3.1 from Greene (1989). Given a matroid M = (E, F) and some subset A of E, finds a minimal spanning subset of A. If A = E, this finds a basis of M. If A is a basis, this finds A.
"""
minimal_spanning_subset(M::ClosedSetsMatroid, A::Integer) = _mss(M, 0, A)
minimal_spanning_subset(M::ClosedSetsMatroid, A) = minimal_spanning_subset(M, set_to_bits(A))

minimal_spanning_subset(M::FullMatroid, A::Integer) = _mss(M, 0, A)
minimal_spanning_subset(M::FullMatroid, A) = minimal_spanning_subset(M, set_to_bits(A))

function _mss(M, j::Integer, Ā::Integer)
  B = [Ā&F for F in M.F[j+1] if Base.count_ones(Ā&F) > j]

  while length(B) == 0
    if j >= Base.count_ones(Ā)-1 return Ā end
    
    j += 1
    B = [Ā&F for F in M.F[j+1] if Base.count_ones(Ā&F) > j]
  end

  _mss(M, j, Ā&~rand_el(reduce(|, B)))
end



