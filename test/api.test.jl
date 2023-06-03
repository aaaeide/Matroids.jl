import LinearAlgebra: I

@testset "Exchange graph" begin
  # Every agent likes every item.
  matroids = [FreeMatroid(5) for _ in 1:5]
  # Every agent has 1 item.
  A = Matrix{Bool}(I, 5, 5) |> BitMatrix
  # Every item should be exchangeable with every other item.
  @test exchange_graph(matroids, A) == complete_digraph(5)
end