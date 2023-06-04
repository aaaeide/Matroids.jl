

@testset "mms_i" begin
  V = MatroidRank([FreeMatroid(3) for _ in 1:3], 3)
  for i in 1:3 @test mms_i(V, i) == 1 end
  
  V = MatroidRank([FreeMatroid(3) for _ in 1:4], 3)
  for i in 1:3 @test mms_i(V, i) == 0 end
end