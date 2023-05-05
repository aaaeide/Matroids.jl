struct ClosedSetsMatroid{T}
  n::Integer # Size of universe
  r::Integer # Final rank (r == length(F)).
  F::Vector{Set{T}} # Closed sets by rank
  rank::Dict{T, UInt8} # Mapping from sets to rank.
  Type::DataType
end

struct FullMatroid{T}
  n::Integer
  r::Integer
  F::Vector{Set{T}} # Closed sets by rank
  I::Vector{Set{T}} # Independent sets by rank
  C::Set{T} # Circuits
  rank::Dict{T, UInt8}
  Type::DataType
end

struct UniformMatroid
  n::Integer
  r::Integer
end

FreeMatroid(n) = UniformMatroid(n, n)


Matroid = Union{ClosedSetsMatroid, FullMatroid, UniformMatroid}