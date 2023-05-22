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

@testset "matroid properties" begin
  include("../src/kmc.jl")
  # The example from Knuth (1974) section 3.
  n = 10
  enlargements = [nothing, [0x1a, 0x222, 0x64, 0x128, 0x288, 0x10c]]
  M = knuth_matroid(n, enlargements)
  
  B = Set([0x00a9,  0x0131,  0x0207,  0x00c3,  0x008d,  0x0123,  0x00c5,  0x0087,  0x02a1,  0x01c1,  0x0185,  0x0291,  0x0099,  0x0113,  0x0151,  0x001d,  0x0213,  0x02c1,  0x010b,  0x008b,  0x020d,  0x0189,  0x0039,  0x00e1,  0x0285,  0x01a1,  0x0311,  0x0191,  0x0245,  0x0381,  0x0107,  0x00a3,  0x00d1,  0x0143,  0x0231,  0x00b1,  0x004b,  0x0243,  0x0017,  0x0035,  0x0053,  0x0115,  0x0261,  0x002b,  0x020b,  0x00c9,  0x0063,  0x0055,  0x0059,  0x0119,  0x0249,  0x0183,  0x0229,  0x0251,  0x0027,  0x0033,  0x0309,  0x0225,  0x0219,  0x0321,  0x0095,  0x0047,  0x0283,  0x0071,  0x000f,  0x0305,  0x0341,  0x0215,  0x0093,  0x00a5,  0x0303])
  
  C = Set([0x001a, 0x0222, 0x002c, 0x004c, 0x010c, 0x0064, 0x0124, 0x0144, 0x0068, 0x0128, 0x0148, 0x0288, 0x0160, 0x008e, 0x020e, 0x0036, 0x0056, 0x0096, 0x0116, 0x0216, 0x00a6, 0x00c6, 0x0246, 0x0186, 0x0286, 0x0306, 0x00aa, 0x00ca, 0x024a, 0x018a, 0x030a, 0x0072, 0x00b2, 0x0132, 0x00d2, 0x0152, 0x0252, 0x0192, 0x0292, 0x0312, 0x00e2, 0x01a2, 0x01c2, 0x02c2, 0x0342, 0x0382, 0x009c, 0x021c, 0x00b4, 0x0234, 0x00d4, 0x0254, 0x0194, 0x0294, 0x0314, 0x02a4, 0x02c4, 0x0384, 0x00b8, 0x0238, 0x00d8, 0x0258, 0x0198, 0x0318, 0x00f0, 0x0270, 0x01b0, 0x02b0, 0x0330, 0x01d0, 0x02d0, 0x0350, 0x0390, 0x02e0, 0x03a0, 0x03c0])

  for base in B @test is_indep(M, base) == true end
  for circ in C @test is_indep(M, circ) == false end

  for base in B @test is_circuit(M, base) == false end
  for circ in C @test is_circuit(M, circ) == true end

  @test minimal_spanning_subset(M, 2^n-1) in B
  @test minimal_spanning_subsets(M, UInt16(2^n-1)) == B
  @test bases(M) == B
end