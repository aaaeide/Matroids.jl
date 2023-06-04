using Matroids
using Test

@testset "Matroids.jl" begin
    include("../src/types.jl")
    include("../src/bitset_utils.jl")

    include("../src/properties.jl")
    include("../src/api.jl")
    
    include("../src/constructions/knuth.jl")
    include("../src/constructions/erect.jl")

    include("../src/fair_allocation/mrf_profile.jl")
    include("../src/fair_allocation/algorithms.jl")
    include("../src/fair_allocation/fairness.jl")

    @test Matroids.greet() == "A matroid /ˈmeɪtrɔɪd/ is a structure that abstracts and generalizes the notion of linear independence in vector spaces."
    
    include("constructions/knuth.test.jl")
    include("constructions/erect.test.jl")

    include("properties.test.jl")
    include("api.test.jl")

    include("bitset_utils.test.jl")

    include("fair_allocation/algorithms.test.jl")
    include("fair_allocation/fairness.test.jl")
end
