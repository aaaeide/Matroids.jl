using Graphs
using Allocations

@testset "Exchange graphs and transfer paths" begin
  # Every agent likes every item.
  matroids = [FreeMatroid(5) for _ in 1:5]
  
  # Every agent has 1 item.
  A = Allocation(5,5)
  for i in 1:5
    give!(A, i, i)
  end

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


  # Three small matroids that require some transfers on the exchange
  # graph to get a partition of the ground set.

  # The ground set:       E = [1       2       3       4       5]
  g1 = SimpleGraph{Int64}(5, [[2, 3], [1, 3], [1, 2], [4],    [5]])
  g2 = SimpleGraph{Int64}(5, [[1],    [3, 4], [2, 4], [2, 3], [5]])
  g3 = SimpleGraph{Int64}(5, [[1],    [2],    [4, 5], [3, 5], [3, 4]])

  ms = [GraphicMatroid(g1), GraphicMatroid(g2), GraphicMatroid(g3)]

  (partition, junk) = matroid_partition_knuth73(ms)
  for (i, set) in enumerate(partition)
    @test is_indep(ms[i], set) || "set $set not indep in ms[$i]"
  end
end
