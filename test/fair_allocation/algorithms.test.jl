using Allocations

@testset "Yankee Swap" begin
  # n = m = 4
  # V = MatroidRank([FreeMatroid(m) for _ in 1:n], m)
  # A = yankee_swap(V)
  # for i in 1:n
  #   @test value(V, i, A) == 1
  # end

  # V = MatroidRank([random_knuth_matroid(16, [0,6]) for _ in 1:7], 16)
  # A = yankee_swap(V)
  # @test check_ef1(V, A)
  # @test check_efx(V, A)
  # @test check_efx0(V, A)
end