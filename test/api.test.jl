using Test
using Graphs

include("../src/api.jl")

@testset "bfs" begin
  # Zachary's karate club
  karate = smallgraph(:karate)
  @test bfs(karate, 1, x -> x==18) == [1,18]
  @test bfs(karate, 1, x->x==34) |> length == 3
  @test bfs(karate, 1, x->x==16) == [1,3,33,16]
end