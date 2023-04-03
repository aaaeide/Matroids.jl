"""
    set_to_bits(set, T=UInt64)

Converts a Set to its bitstring representation.
"""
function set_to_bits(set, T=UInt64)
  T(sum(2^x for x in set))
end

# Drawing graphs
using GraphPlot
using Compose
import Fontconfig, Cairo

function draw_graph(g, nodelabel, filename="graph.png")
  draw(PNG(filename, 32cm, 32cm), gplot(g, nodelabel=nodelabel))
end