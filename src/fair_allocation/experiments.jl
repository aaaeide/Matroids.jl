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
function experiments(gen_profile, algo, n)
  results = Dict("α-EF1" => [],
                 "α-EFX_+" => [],
                 "α-EFX_0" => [],
                 "α-PROP" => [],
                 "α-PROPX_+" => [],
                #  "α-PROPX_0" => [],
                 "α-MMS" => [],
                 "time" => [],
                 "bytes" => [])

  println(" Time      | Bytes allocated | α-EF1 | α-EFX_+ | α-EFX_0 | α-PROP | α-PROPX_+ | α-PROPX_0 | α-MMS ")
  println("-----------|-----------------|-------|---------|---------|--------|----------|-----------|-------")
  
  for _ in 1:n
    V = gen_profile()
    res = @timed A = algo(V)

    push!(results["α-EF1"], check_alpha_ef_(V, A, value_1))
    push!(results["α-EFX_+"], check_alpha_ef_(V, A, value_x))
    push!(results["α-EFX_0"], check_alpha_ef_(V, A, value_x0))
    push!(results["α-PROP"], check_alpha_thresh(V, A, prop))
    push!(results["α-PROPX_+"], check_alpha_thresh(V, A, prop_x))
    # push!(results["α-PROPX_0"], check_alpha_thresh(V, A, prop_x0))
    push!(results["α-MMS"], check_alpha_mms(V, A))
    push!(results["time"], res.time)
    push!(results["bytes"], res.bytes)

    for (k, v) in results println("$k\t\t$(v[end])") end
    println("\n")
  end
end

"""
Generate a MatroidRank valuation profile with n matroids over m elements. 

The matroids are generated using random_knuth_matroid(n,p,T).
"""
function gen_closedsets_matroidrank_profile(n, m, p, T)
  function gen()
    ms = Array{ClosedSetsMatroid}(undef, n)

    Threads.@threads for i in 1:n
      ms[i] = random_knuth_matroid(m, p, T)
    end

    return MatroidRank(ms, m)
  end

  return gen
end

