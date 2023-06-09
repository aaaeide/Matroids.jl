using Allocations

"""
Given a function to generate a random valuation profile and an algorithm to test, runs the algorithm on n generated valuation profiles, checking the fairness of each allocation.

Fairness notions checked:
- EF1
- EFX_+
- EFX_0
- PROP
- PROPX_+
- PROPX_0
- MMS

Also: average time to compute and average bytes allocated.
"""
function experiments(profile_generators, algo, n)
 
  for (name, gen_profile) in profile_generators
    results = Dict{String, Union{Vector{Any}, Float64}}("α-EF1" => [],
                   "α-EFX_+" => [],
                   "α-EFX_0" => [],
                   "α-PROP" => [],
                   "α-PROPX_+" => [],
                   "α-PROPX_0" => [],
                   "α-MMS" => [],
                   "time_alloc" => [],
                   "bytes_alloc" => [],
                   "time_gener" => [],
                   "bytes_gener" => [])

    println("\n\n$name")

    for _ in 1:n
      res_gener = @timed V = gen_profile()
      res_alloc = @timed A = algo(V)
  
      push!(results["α-EF1"], 
          check_alpha_ef_(V, A, value_1; indep=true))

      push!(results["α-EFX_+"], 
          check_alpha_ef_(V, A, value_x; indep=true))

      push!(results["α-EFX_0"], 
          check_alpha_ef_(V, A, value_x0; indep=true))

      push!(results["α-PROP"], 
          check_alpha_thresh(V, A, prop; indep=true))

      push!(results["α-PROPX_+"], 
          check_alpha_thresh(V, A, prop_x; indep=true))

      push!(results["α-PROPX_0"], 
          check_alpha_thresh(V, A, prop_x0; indep=true))

      push!(results["α-MMS"], 
          check_alpha_mms(V, A; indep=true))

      push!(results["time_alloc"], res_alloc.time)
      push!(results["bytes_alloc"], res_alloc.bytes)
      
      push!(results["time_gener"], res_gener.time)
      push!(results["bytes_gener"], res_gener.bytes)
    end
    
    for (k,v) in results
      results[k] = sum(v)/length(v)
  
      if k ∉ ["time_alloc", "time_gener"]
        results[k] = trunc(results[k], digits=3)
      end
    end
    
    time_alloc = format_time(results["time_alloc"])
    time_gener = format_time(results["time_gener"])

    bytes_alloc = Base.format_bytes(results["bytes_alloc"])
    bytes_gener = Base.format_bytes(results["bytes_gener"])
        
    println("Time (Gen.) | Bytes (Gen.) | Time (Alloc.) | Bytes (Alloc.) | α-EF1 | α-EFX_+ | α-EFX_0 | α-PROP | α-PROPX_+ | α-PROPX_0 | α-MMS ")
    println("------------|--------------|---------------|----------------|-------|---------|---------|--------|-----------|-----------|-------")
    println("$(rpad(time_gener, 11, " ")) | $(rpad(bytes_gener, 12, " ")) | $(rpad(time_alloc, 13, " ")) | $(rpad(bytes_alloc, 14, " ")) | $(rpad(results["α-EF1"], 5, " ")) | $(rpad(results["α-EFX_+"], 7, " ")) | $(rpad(results["α-EFX_0"], 7, " ")) | $(rpad(results["α-PROP"], 6, " ")) | $(rpad(results["α-PROPX_+"], 9, " ")) | $(rpad(results["α-PROPX_0"], 9, " ")) | $(results["α-MMS"])")
  end
end

function format_time(secs)
  time = secs * 1e9
  if time > 1e9
    time = "$(round(time / 1e9, digits=3))s"
  elseif time > 1e6
    time = "$(round(time / 1e6, digits=3))ms"
  elseif time > 1e3
    time = "$(round(time / 1e3, digits=3))µs"
  else
    time = "$(round(time, digits=3))ns"
  end

  return time
end


profile_gens = [
  ("8 x GraphicMatroid(random_er_graph(64))", 
  gen_matroidrank_profile(8, 
    () -> GraphicMatroid(random_er_graph(64)), 
    GraphicMatroid)),
  ("8 x GraphicMatroid(random_ws_graph(64))", 
  gen_matroidrank_profile(8, 
    () -> GraphicMatroid(random_ws_graph(64)), 
    GraphicMatroid)),
  ("8 x GraphicMatroid(random_ba_graph(64))", 
  gen_matroidrank_profile(8, 
    () -> GraphicMatroid(random_ba_graph(64)), 
    GraphicMatroid)),
  ("4 x random_knuth_matroid(13, [0,6,0])", 
  gen_matroidrank_profile(4, 
    () -> random_knuth_matroid(13, [0,6,0]), 
    ClosedSetsMatroid)),
  ("4 x random_knuth_matroid(16, [0,7,2])", 
  gen_matroidrank_profile(4, 
    () -> random_knuth_matroid(16, [0,7,2]), 
    ClosedSetsMatroid)),
  ("4 x random_knuth_matroid(19, [0,8,3,1])", 
  gen_matroidrank_profile(4, 
    () -> random_knuth_matroid(19, [0,8,3,1], UInt32), 
    ClosedSetsMatroid)),
]