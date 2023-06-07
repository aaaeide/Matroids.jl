
@testset "ef1" begin
  V = MatroidRank([UniformMatroid(10, 6), UniformMatroid(10, 6)], 10)
  @test value_1(V, 1, 1:7) == 6
end

@testset "mms_i" begin
  V = MatroidRank([FreeMatroid(3) for _ in 1:3], 3)
  for i in 1:3 @test mms_i(V, i) == 1 end
  
  V = MatroidRank([FreeMatroid(3) for _ in 1:4], 3)
  for i in 1:3 @test mms_i(V, i) == 0 end
end

@testset "alpha-ef_" begin
  V = MatroidRank([UniformMatroid(5, 2), UniformMatroid(5, 2), UniformMatroid(5, 2)], 5)
  A = alloc_algmms_bv21(V)
  @test check_alpha_ef_(V, A, value) == 0.5
  @test check_alpha_ef_(V, A, value_1) == 1
end