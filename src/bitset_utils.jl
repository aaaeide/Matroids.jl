"""
    set_to_bits(set, T=UInt64)

Converts a Set to its bitstring representation. 

NB! set_to_bits and its inverse bits_to_set assume the elements (in the set representation) start at 1.
"""
function set_to_bits(set, T=UInt64)
  if length(set) == 0 return T(0) end
  T(reduce(+, (2^(x-1) for x in set), init=0))
end


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


function rand_el(S::Integer)
  x = rand([2^(i-1) for (i,c) in enumerate(reverse(bitstring(S))) if c == '1'])
  return convert(typeof(S), x)
end


"""
    add_el(S::Integer, e::Integer)

Returns the set S ∪ {e}. NB! If S is an n-digit set, this assumes 0≤e≤n-1.
"""
add_el(S::Integer, e::Integer) = S | 1<<e

# # Drawing graphs
# using GraphPlot
# using Compose
# import Fontconfig, Cairo

# function draw_graph(g, nodelabel, filename="graph.png")
#   draw(PNG(filename, 64cm, 64cm), gplot(g, nodelabel=nodelabel))
# end