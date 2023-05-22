using Graphs

struct GraphicMatroid
  g::Graph
  n::Integer
  r::Integer
  GraphicMatroid(g::Graph) = new(g, ne(g), length(kruskal_mst(g)))
end

function rank(m::GraphicMatroid, S)
  edgelist = [e for (i, e) in enumerate(edges(m.g)) if i in S]
  subgraph, _vmap = induced_subgraph(m.g, edgelist)
  return length(kruskal_mst(subgraph))
end

function is_indep(m::GraphicMatroid, S)
  edgelist = [e for (i, e) in enumerate(edges(g)) if i in S]
  subgraph, _vmap = induced_subgraph(m.g, edgelist)
  return !is_cyclic(subgraph)
end

function is_circuit(m::GraphicMatroid, S)
  #TODO
end

function closure(m::GraphicMatroid, S)
  edgelist = [e for (i, e) in enumerate(edges(m.g)) if i in S]
  _sg, vmap = induced_subgraph(m.g, edgelist)
  return [e for e in edges(m.g) if [e.src, e.dst] ⊆ vmap]
end