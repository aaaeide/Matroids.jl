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

    # Builds a family of larger sets
    @test family([[1,2], [3,4]]) == Set([Set([1,2]), Set([3,4])])
  end

  @testset "generate_covers" begin
    elements = [0,1,2,3,4,5,6,7,8,9]
    ground_set = Set(elements)

    @test generate_covers(family([]), ground_set) == family(elements)
    
    # The covers of the singleton subsets is the family of all pairs
    @test generate_covers(family(elements), ground_set) == 
      Set([Set([i,j]) for i in ground_set for j in ground_set if i != j])
  end 

  @testset "enlarge" begin
    @test enlarge(family([1,2,3]), [Set(e) for e in [7,8,9]]) == family([1,2,3,7,8,9])
  end

  @testset "superpose" begin
    @test superpose(family([1,2,3]), family([])) == family([1,2,3])

    # This sets up the example from Knuth (1974) p. 345
    elements = [0,1,2,3,4,5,6,7,8,9]
    ground_set = Set(elements)
    F1 = family(elements)
    F2 = generate_covers(family(elements), ground_set)
    F2 = enlarge(F2, family([
      [1,3,4], [1,5,9], [2,5,6], [3,5,8], [3,7,9], [2,3,8]]
    ))

    @test superpose(F2, F1) == family([
      [0,1], [0,2], [0,3], [0,4], [0,5], [0,6], [0,7], [0,8], [0,9], 
      [1,2], [1,3,4], [1,5,9], [1,6], [1,7], [1,8], 
      [2,3,5,6,8], [2,4], [2,7], [2,9],
      [3,7,9],
      [4,5], [4,6], [4,7], [4,8], [4,9],
      [5,7],
      [6,7], [6,9],
      [7,8],
      [8,9]])
  end
end