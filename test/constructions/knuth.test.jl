using Test

@testset "Pi-based KMC example" begin
  # The example from Knuth (1974) section 3.
  n = 10
  enlargements = [nothing, [0x1a, 0x222, 0x64, 0x128, 0x288, 0x10c]]
  
  F0 = Set(0) # r=0: The empty set alone.
  F1 = Set([set_to_bits(Set(i)) for i in 1:n]) # r=1: Singleton subsets of E.
  F2 = Set([set_to_bits(set) for set in 
    [[1, 2], [1, 3], [1, 4], [1, 5], [1, 6], [1, 7], [1, 8], [1, 9], 
    [1, 10], [2, 3], [2, 4, 5], [2, 6, 10], [2, 7], [2, 8], [2, 9], 
    [3, 4, 6, 7, 9], [3, 5], [3, 8], [3, 10], [4, 8, 10], [5, 6], 
    [5, 7], [5, 8], [5, 9], [5, 10], [6, 8], [7, 8], [7, 10], 
    [8, 9], [9, 10]]])
  F3 = Set([set_to_bits(set) for set in 
    [[1, 2, 3], [1, 2, 4, 5], [1, 2, 6, 10], [1, 2, 7], [1, 2, 8], 
    [1, 2, 9], [1, 3, 4, 6, 7, 9], [1, 3, 5], [1, 3, 8], [1, 3, 10], 
    [1, 4, 8, 10], [1, 5, 6], [1, 5, 7], [1, 5, 8], [1, 5, 9], 
    [1, 5, 10], [1, 6, 8], [1, 7, 8], [1, 7, 10], [1, 8, 9], 
    [1, 9, 10], [2, 3, 4, 5, 6, 7, 8, 9, 10]]])
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