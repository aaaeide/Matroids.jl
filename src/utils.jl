"""
    set_to_bits(set, T=UInt64)

Converts a Set to its bitstring representation.
"""
function set_to_bits(set::Set, T=UInt64)
  if length(set) == 0 return T(0) end
  T(sum(2^x for x in set))
end

set_to_bits(vec, T=UInt64) = set_to_bits(Set(vec), T)

"""
    function bits_to_set(bits)

Converts a set represented as a bitstring to a Set.
"""
function bits_to_set(bits)
  Set(i-1 for (i, c) in enumerate(reverse(bitstring(bits))) if c == '1')
end



# # Drawing graphs
# using GraphPlot
# using Compose
# import Fontconfig, Cairo

# function draw_graph(g, nodelabel, filename="graph.png")
#   draw(PNG(filename, 64cm, 64cm), gplot(g, nodelabel=nodelabel))
# end