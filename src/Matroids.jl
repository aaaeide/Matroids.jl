module Matroids

# Write your package code here.

export greet

greet() = "A matroid /ˈmeɪtrɔɪd/ is a structure that abstracts and generalizes the notion of linear independence in vector spaces."

include("types.jl")
include("api.jl")
include("bitset_utils.jl")
include("axioms.jl")
include("properties.jl")

include("constructions/knuth.jl")
include("constructions/erect.jl")

include("fair_allocation/mrf_profile.jl")
include("fair_allocation/fairness.jl")
include("fair_allocation/algorithms.jl")


export Matroid, ZeroMatroid, FreeMatroid, UniformMatroid, GraphicMatroid, MatroidRank, yankee_swap
end
