using Test

include("../src/properties.jl")
include("../src/kmc.jl")

@testset "KMC v1: Bases, rank and closure" begin
  # The example from Knuth (1974) section 3
  E = Set([0,1,2,3,4,5,6,7,8,9])
  enlargements = [nothing, family([
    [1,3,4], [1,5,9], [2,5,6], [3,5,8], [3,7,9], [2,3,8]]
  )]

  M = knuth_matroid_construction_v1(E, enlargements)

  @test_throws ArgumentError closure(M, [-1, -2])

  # In the Knuth example, the singleton sets are all closed.
  @test rank_and_closure(M, [0]) == (1, Set([0]))

  # Closed sets of rank 2
  @test rank_and_closure(M, [1,2]) == (2, Set([1,2]))
  @test rank_and_closure(M, [1,5]) == (2, Set([1,5,9]))
  @test rank_and_closure(M, [2,3]) == (2, Set([2,3,5,6,8]))

  # Closed sets of rank 3
  @test rank_and_closure(M, [0,1,5]) == (3, Set([0,1,5,9]))
  @test rank_and_closure(M, [1,2,3]) == (3, Set([1,2,3,4,5,6,7,8,9]))

  # Closed sets of rank 4 (the example matroid is rank 4, so this is E).
  @test rank_and_closure(M, [0,1,2,3]) == (4, E)

  @test bases(M) == family([
      [0,1,2,3], [0,1,2,4], [0,1,2,5], [0,1,2,6], [0,1,2,7], [0,1,2,8], [0,1,2,9],
      [0,1,3,5], [0,1,3,6], [0,1,3,7], [0,1,3,8], [0,1,3,9],
      [0,1,4,5], [0,1,4,6], [0,1,4,7], [0,1,4,8], [0,1,4,9],
      [0,1,5,6], [0,1,5,7], [0,1,5,8],
      [0,1,6,7], [0,1,6,8], [0,1,6,9],
      [0,1,7,8], [0,1,7,9],
      [0,1,8,9],
      [0,2,3,4], [0,2,3,7], [0,2,3,9],
      [0,2,4,5], [0,2,4,6], [0,2,4,7], [0,2,4,8], [0,2,4,9],
      [0,2,5,7], [0,2,5,9],
      [0,2,6,7], [0,2,6,9],
      [0,2,7,8], [0,2,7,9],
      [0,2,8,9], 
      [0,3,4,5], [0,3,4,6], [0,3,4,7], [0,3,4,8], [0,3,4,9],
      [0,3,5,7], [0,3,5,9],
      [0,3,6,7], [0,3,6,9],
      [0,3,7,8],
      [0,3,8,9],
      [0,4,5,6], [0,4,5,7], [0,4,5,8], [0,4,5,9],
      [0,4,6,7], [0,4,6,8], [0,4,6,9],
      [0,4,7,8], [0,4,7,9],
      [0,4,8,9],
      [0,5,6,7], [0,5,6,9],
      [0,5,7,8], [0,5,7,9],
      [0,5,8,9],
      [0,6,7,8], [0,6,7,9], 
      [0,6,8,9],
      [0,7,8,9],
    ])
end