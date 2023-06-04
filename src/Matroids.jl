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

export Matroid, ZeroMatroid, FreeMatroid, UniformMatroid, ClosedSetsMatroid, GraphicMatroid
export is_indep, rank, is_circuit, minimal_spanning_subset, minimal_spanning_subsets, bases, closure

export knuth_matroid, random_knuth_matroid, erect_v1

export MatroidRank, na, ni, value
export check_ef, check_ef1, check_efx, check_efx0, check_prop, check_prop1, check_propx, check_propx0, mms_i
export alloc_yankee_swap_vz22, alloc_bciz21
end
