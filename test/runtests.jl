using Matroids
using Test

@testset "Matroids.jl" begin
    include("../src/utils.jl")
    include("../src/types.jl")
    include("../src/kmc.jl")
    include("../src/erect.jl")
    include("../src/properties.jl")

    # Write your tests here.

    @test Matroids.greet() == "A matroid /ˈmeɪtrɔɪd/ is a structure that abstracts and generalizes the notion of linear independence in vector spaces."

    
    include("kmc.test.jl")
    include("erect.test.jl")
    include("properties.test.jl")
end
