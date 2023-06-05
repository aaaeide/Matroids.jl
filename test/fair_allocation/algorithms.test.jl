using Allocations

@testset "bv_to_bundle" begin
  @test bv_to_bundle(BitVector([1,1,0,0,1])) == Set([1,2,5])
end

@testset "Algorithms" begin
  n = m = 4
  V1 = MatroidRank([FreeMatroid(m) for _ in 1:n], m)
  V2 = MatroidRank([random_knuth_matroid(16, [0,6,3,2]) for _ in 1:7], 16)

  # This scenario results in some non-trivial transfer paths, and includes 
  # an element no agent likes.
  g1 = SimpleGraph{Int64}(6, [[2, 3], [1, 3], [1, 2], [4],    [5],    [6]])
  g2 = SimpleGraph{Int64}(6, [[1],    [3, 4], [2, 4], [2, 3], [5],    [6]])
  g3 = SimpleGraph{Int64}(6, [[1],    [2],    [4, 5], [3, 5], [3, 4], [6]])

  ms = [GraphicMatroid(g1), GraphicMatroid(g2), GraphicMatroid(g3)]

  V3 = MatroidRank(ms, 6)
  
  @testset "Yankee Swap" begin
    A = alloc_yankee_swap_vz22(V1)
    for i in 1:n
      @test value(V1, i, A) == 1
    end
    
    A = alloc_yankee_swap_vz22(V2)
    @test check_ef1(V2, A)
    @test check_efx(V2, A)
    @test check_efx0(V2, A)

    A = alloc_yankee_swap_vz22(V3)
    @test check_ef1(V3, A)
    @test check_efx(V3, A)
    @test check_efx0(V3, A)
  end
  


  @testset "Envy-induced transfer" begin
    A = alloc_eit_bciz21(V1)
    for i in 1:n
      @test value(V1, i, A) == 1
    end
    
    A = alloc_eit_bciz21(V2)
    @test check_ef1(V2, A)
    @test check_efx(V2, A)
    @test check_efx0(V2, A)
    
    A = alloc_eit_bciz21(V3)
    @test check_ef1(V3, A)
    @test check_efx(V3, A)
    @test check_efx0(V3, A)
  end
  
  
  
  @testset "AlgMMS" begin
    A = alloc_algmms_bv21(V1)
    for i in 1:n
      @test value(V1, i, A) == 1
    end

    A = alloc_algmms_bv21(V2)
    @test check_mms(V2, A)
    
    A = alloc_algmms_bv21(V3)
    @test check_mms(V3, A)
  end
end
