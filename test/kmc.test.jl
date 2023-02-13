using Test
include("../src/kmc.jl")

@testset "Set-based KMC" begin
  @testset "family helper" begin
    # Builds the family consisting of the empty set alone
    @test family([]) == Set([Set([])])

    # Builds the family of all singleton subsets
    @test family([1,2,3]) == Set([Set(e) for e in [1,2,3]])

    # Builds a family of larger sets
    @test family([[1,2], [3,4]]) == Set([Set([1,2]), Set([3,4])])
  end

  @testset "generate_covers_v1" begin
    elements = [0,1,2,3,4,5,6,7,8,9]
    ground_set = Set(elements)

    @test generate_covers_v1(family([]), ground_set) == family(elements)
    
    # The covers of the singleton subsets is the family of all pairs
    @test generate_covers_v1(family(elements), ground_set) == 
      Set([Set([i,j]) for i in ground_set for j in ground_set if i != j])
  end

  @testset "superpose_v1" begin
    @test superpose_v1!(family([1,2,3]), family([])) == family([1,2,3])

    # This sets up (part of) the example from Knuth (1974) p. 345
    elements = [0,1,2,3,4,5,6,7,8,9]
    ground_set = Set(elements)
    F1 = family(elements)
    F2 = generate_covers_v1(family(elements), ground_set) ∪ family([
      [1,3,4], [1,5,9], [2,5,6], [3,5,8], [3,7,9], [2,3,8]
    ])

    @test superpose_v1!(F2, F1) == family([
      [0,1], [0,2], [0,3], [0,4], [0,5], [0,6], [0,7], [0,8], [0,9], 
      [1,2], [1,3,4], [1,5,9], [1,6], [1,7], [1,8], 
      [2,3,5,6,8], [2,4], [2,7], [2,9],
      [3,7,9],
      [4,5], [4,6], [4,7], [4,8], [4,9],
      [5,7],
      [6,7], [6,9],
      [7,8],
      [8,9]
    ])
  end

  @testset "knuth_matroid_construction_v1" begin
    # The example from Knuth (1974) section 3
    E = Set([0,1,2,3,4,5,6,7,8,9])
    enlargements = [nothing, family([
      [1,3,4], [1,5,9], [2,5,6], [3,5,8], [3,7,9], [2,3,8]]
    )]

    # The closed sets, grouped by rank
    F0 = family([])
    F1 = family(E) # The family of singleton sets
    F2 = family([
      [0,1], [0,2], [0,3], [0,4], [0,5], [0,6], [0,7], [0,8], [0,9], 
      [1,2], [1,3,4], [1,5,9], [1,6], [1,7], [1,8], 
      [2,3,5,6,8], [2,4], [2,7], [2,9],
      [3,7,9],
      [4,5], [4,6], [4,7], [4,8], [4,9],
      [5,7],
      [6,7], [6,9],
      [7,8],
      [8,9]
    ])
    F3 = family([
      [0,1,2], [0,1,3,4], [0,1,5,9], [0,1,6], [0,1,7], [0,1,8],
      [0,2,3,5,6,8], [0,2,4], [0,2,7], [0,2,9],
      [0,3,7,9], 
      [0,4,5], [0,4,6], [0,4,7], [0,4,8], [0,4,9],
      [0,5,7],
      [0,6,7], [0,6,9],
      [0,7,8],
      [0,8,9],
      [1,2,3,4,5,6,7,8,9]
    ])
    F4 = family([E]) # The family of only one set - E
    F = [F0, F1, F2, F3, F4]

    @test knuth_matroid_construction_v1(E, enlargements) == (E, F)
  end
end

@testset "Bitwise KMC" begin
  @testset "generate_covers_v2" begin
    # The covers of the empty set are the singleton sets.
    @test generate_covers_v2(Set(0), 4) == Set([0x01, 0x02, 0x04, 0x08]) 
    
    # The covers of the singleton sets is the family of all pairs.
    @test generate_covers_v2(Set([1 << i for i in 0:3]), 4) == 
      Set([(1 << i) | (1 << j) for i in 0:3 for j in 0:3 if i != j])
  end

  @testset "bitwise_superpose!" begin
    # This sets up (part of) the example from Knuth (1974) p. 345
    n = 10
    F0 = Set(0) # r=0: The empty set alone.
    F1 = Set([1 << i for i in 0:n-1]) # r=1: Singleton subsets of E.
    F2 = generate_covers_v2(F1, n) ∪ [set_to_bits(set) for set in [[1,3,4], [1,5,9], [2,5,6], [3,5,8], [3,7,9], [2,3,8]]]
    
    @test bitwise_superpose!(F1, F0) == F1

    result = bitwise_superpose!(F2, F1)
    expect = Set([set_to_bits(set) for set in [
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

    @test result == expect
    @test setdiff(result, expect) == Set()
    @test setdiff(expect, result) == Set()
  end

  @testset "knuth_matroid_construction_v2" begin
    # The example form Knuth (1974) section 3.
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

    @test knuth_matroid_construction_v2(n, enlargements) == (n, F)
  end
end

@testset "Sorted bitwise KMC" begin
  @testset "sorted_bitwise_superpose!" begin
    # This sets up (part of) the example from Knuth (1974) p. 345
    n = 10
    F0 = Set(0) # r=0: The empty set alone.
    F1 = Set([1 << i for i in 0:n-1]) # r=1: Singleton subsets of E.
    F2 = generate_covers_v2(F1, n) ∪ [set_to_bits(set) for set in [[1,3,4], [1,5,9], [2,5,6], [3,5,8], [3,7,9], [2,3,8]]]
    
    @test sorted_bitwise_superpose!(F1, F0) == F1

    result = sorted_bitwise_superpose!(F2, F1)
    expect = Set([set_to_bits(set) for set in [
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

    @test result == expect
    @test setdiff(result, expect) == Set()
    @test setdiff(expect, result) == Set()
  end

  @testset "knuth_matroid_construction_v3" begin
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

    @test knuth_matroid_construction_v3(n, enlargements) == (n, F)
  end
end