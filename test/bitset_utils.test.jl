@testset "bitset utils" begin
  @test set_to_bits([1,3,4]) == 0b00001101
  @test bits_to_set(0x16c) == Set([3,4,6,7,9])

  @test [3,4,6,7,9] |> set_to_bits |> bits_to_set == Set([3,4,6,7,9])
  @test Set([3,4,6,7,9]) |> set_to_bits |> bits_to_set == Set([3,4,6,7,9])

  @test add_el(0b001, 2) == 0b101
end