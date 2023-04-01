"""
    set_to_bits(set, T=UInt64)

Converts a Set to its bitstring representation.
"""
function set_to_bits(set, T=UInt64)
  T(sum(2^x for x in set))
end
