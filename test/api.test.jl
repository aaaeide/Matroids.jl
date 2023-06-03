using Graphs
import LinearAlgebra: I

@testset "Exchange graphs and transfer paths" begin
  # Every agent likes every item.
  matroids = [FreeMatroid(5) for _ in 1:5]
  # Every agent has 1 item.
  A = Matrix{Bool}(I, 5, 5) |> BitMatrix
  # Every item should be exchangeable with every other item.
  @test exchange_graph(matroids, A) == complete_digraph(5)

  # Zachary's Karate Club.
  k = smallgraph(:karate)
  @test find_shortest_path(k, [1], [16]) |> length == 4
  @test find_shortest_path(k, [1,2,3], [1,18,4]) == [1]

  # A graph with no edges has no transfer paths.
  g = SimpleGraph(10, 0)
  @test find_shortest_path(g, [1,2,3], [4,5,6]) === nothing
  @test find_shortest_path(g, [1,2,3,4], [4,5]) == [4]
end

@testset "Transfer operation" begin

  # A = 1 0 0 0
  #     0 1 0 0
  #     0 0 1 0
  #     0 0 0 1

  A = Matrix{Bool}(I,4,4) |> BitMatrix
  P = [2,3,4]
  i = 1

  @test transfer!(A, i, P) == [1 1 0 0;
                               0 0 1 0; 
                               0 0 0 1;
                               0 0 0 0]
end