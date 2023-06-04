using Allocations

@testset "bv_to_bundle" begin
  @test bv_to_bundle(BitVector([1,1,0,0,1])) == Set([1,2,5])
end

@testset "Yankee Swap" begin
  n = m = 4
  V = MatroidRank([FreeMatroid(m) for _ in 1:n], m)
  A = alloc_yankee_swap_vz22(V)

  for i in 1:n
    @test value(V, i, A) == 1
  end

  V = MatroidRank([random_knuth_matroid(16, [0,6,3,2]) for _ in 1:7], 16)
  A = alloc_yankee_swap_vz22(V)
  @test check_ef1(V, A)
  @test check_efx(V, A)
  @test check_efx0(V, A)
end

@testset "Envy-induced transfer" begin
  n = m = 4
  V = MatroidRank([FreeMatroid(m) for _ in 1:n], m)
  A = alloc_eit_bciz21(V)

  for i in 1:n
    @test value(V, i, A) == 1
  end
end