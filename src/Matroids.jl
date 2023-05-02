module Matroids

# Write your package code here.

export greet

greet() = "A matroid /ˈmeɪtrɔɪd/ is a structure that abstracts and generalizes the notion of linear independence in vector spaces."

include("types.jl")
include("kmc.jl")
include("api.jl")
include("utils.jl")
include("axioms.jl")

end
