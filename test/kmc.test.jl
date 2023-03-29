using Test
include("../src/kmc.jl")

function set_to_bits(set)
  sum(2^x for x in set)
end

@testset "Pi-based KMC example" begin
  # The example from Knuth (1974) section 3.
  n = 10
  enlargements = [nothing, [0x1a, 0x222, 0x64, 0x128, 0x288, 0x10c]]
  
  F0 = Set(0) # r=0: The empty set alone.
  F1 = Set([1 << i for i in 0:n-1]) # r=1: Singleton subsets of E.
  F2 = Set([set_to_bits(set) for set in [
    [0,1], [0,2], [0,3], [0,4], [0,5], [0,6], [0,7], [0,8], [0,9], 
    [1,2], [1,3,4], [1,5,9], [1,6], [1,7], [1,8], 
    [2,3,5,6,8], [2,4], [2,7], [2,9],
    [3,7,9],
    [4,5], [4,6], [4,7], [4,8], [4,9],
    [5,7],
    [6,7], [6,9],
    [7,8],
    [8,9]
    ]])
  F3 = Set([set_to_bits(set) for set in [
    [0,1,2], [0,1,3,4], [0,1,5,9], [0,1,6], [0,1,7], [0,1,8],
    [0,2,3,5,6,8], [0,2,4], [0,2,7], [0,2,9],
    [0,3,7,9], 
    [0,4,5], [0,4,6], [0,4,7], [0,4,8], [0,4,9],
    [0,5,7],
    [0,6,7], [0,6,9],
    [0,7,8],
    [0,8,9],
    [1,2,3,4,5,6,7,8,9]
    ]])
  F4 = Set(2^n-1) # r=4: The family of only one set - E
  F = [F0, F1, F2, F3, F4]

  M = knuth_matroid(n, enlargements)
  result = M.F
  
  @test result[1] == F0
  @test result[2] == F1
  @test result[3] == F2
  @test result[4] == F3
  @test result[5] == F4
end