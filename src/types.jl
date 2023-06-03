using Graphs
using Multigraphs

abstract type Matroid end

struct ClosedSetsMatroid{T} <: Matroid
  n::Integer # Size of universe
  r::Integer # Final rank (r == length(F)).
  F::Vector{Set{T}} # Closed sets by rank
  rank::Dict{T, Integer} # Mapping from sets to rank.
  Type::DataType
end

struct FullMatroid{T} <: Matroid
  n::Integer
  r::Integer
  F::Vector{Set{T}} # Closed sets by rank
  I::Vector{Set{T}} # Independent sets by rank
  C::Set{T} # Circuits
  rank::Dict{T, Integer}
  Type::DataType
end

struct UniformMatroid <: Matroid
  n::Integer
  r::Integer
end

struct GraphicMatroid <: Matroid
  g::Graph
  n::Integer
  r::Integer
  GraphicMatroid(g::Graph) = new(g, ne(g), length(kruskal_mst(g)))
end

FreeMatroid(n) = UniformMatroid(n, n)
ZeroMatroid(n) = UniformMatroid(n, 0)
