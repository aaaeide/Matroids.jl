"""
    set_to_bits(set, T=UInt64)

Converts a Set to its bitstring representation. 

NB! set_to_bits and its inverse bits_to_set assume the elements start at 1.  
"""
function set_to_bits(set::Set, T=UInt64)
  if length(set) == 0 return T(0) end
  T(sum(2^(x-1) for x in set))
end

set_to_bits(vec, T=UInt64) = set_to_bits(Set(vec), T)


"""
    function bits_to_set(bits)

Converts a set represented as a bitstring to a Set.
"""
function bits_to_set(bits)
  Set(i for (i, c) in enumerate(reverse(bitstring(bits))) if c == '1')
end

"""
    function slimmest_integer(n)

Get the UInt least bitwidth larger than n. 
"""
function slimmest_uint(n)
  if n <= 8 return UInt8 end
  if n <= 16 return UInt16 end
  if n <= 32 return UInt32 end
  if n <= 64 return UInt64 end
  return UInt128
end

# # Drawing graphs
# using GraphPlot
# using Compose
# import Fontconfig, Cairo

# function draw_graph(g, nodelabel, filename="graph.png")
#   draw(PNG(filename, 64cm, 64cm), gplot(g, nodelabel=nodelabel))
# end