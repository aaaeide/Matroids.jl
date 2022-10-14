using Test
include("../src/knuth74.jl")

function family(set_contents::Vector{Any})::Family{Any}
  Set([Set([])])
end

function family(set_contents::Vector{<:Any})::Family{<:Any}
  Set([Set(e) for e in set_contents])
end

@testset "Knuth's matroid construction" begin
  @testset "family helper" begin
    # Builds the family consisting of the empty set alone
    @test family([]) == Set([Set([])])

    # Builds the family of all singleton subsets
    @test family([1,2,3]) == Set([Set(e) for e in [1,2,3]])
  end

  @testset "generate_covers" begin
    elements = [0,1,2,3,4,5,6,7,8,9]
    ground_set = Set(elements)

    @test generate_covers(family([]), ground_set) == family(elements)
    
    @test generate_covers(family(elements), ground_set) == 
      Set([Set([i,j]) for i in ground_set for j in ground_set if i != j])
  end 

  @testset "enlarge" begin
    @test enlarge(family([1,2,3]), [Set(e) for e in [7,8,9]]) == family([1,2,3,7,8,9])
  end

  @testset "superpose" begin
    @test superpose(family([1,2,3]), family([])) == family([1,2,3])

    
  end
end