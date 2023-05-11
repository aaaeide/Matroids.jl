using BenchmarkTools

function avg(arr)
  return sum(arr) / length(arr)
end

function experiment(kmc, test)
  times = []
  bytes = []
  gctimes = []
  ranks = []
  for _ in 1:test.num
    res = @timed kmc(test.n, test.p, test.T)
    push!(times, res.time)
    push!(bytes, res.bytes)
    push!(gctimes, res.gctime)
    push!(ranks, res.value.r)
  end

  t = round(avg(times), digits=7)
  b = floor(Int, avg(bytes))
  g = round(avg(gctimes), digits=7)
  r = round(avg(ranks), digits=4)

  return (t,b,g,r)
end

function experiments(kmc)
  tests = [
    (n=8,  p=[0,1,0],   num=100, T=UInt16),
    (n=10, p=[0,6,0],   num=100, T=UInt16),
    (n=10, p=[0,4,2],   num=100, T=UInt16),
    (n=10, p=[0,0,6],   num=100, T=UInt16),
    (n=10, p=[0,1,1,1], num=100, T=UInt16),
    (n=13, p=[0,8,0],   num=75,  T=UInt16),
    (n=13, p=[0,5,3],   num=75,  T=UInt16),
    (n=13, p=[0,0,8],   num=75,  T=UInt16),
    (n=16, p=[0,10,0],  num=50,  T=UInt16),
    (n=16, p=[0,6,4],   num=50,  T=UInt16),
    (n=16, p=[0,0,10],  num=50,  T=UInt16),
    (n=24, p=[0,14,0],  num=25,  T=UInt32),
    (n=24, p=[0,11,3],  num=25,  T=UInt32),
    (n=24, p=[0,0,14],  num=25,  T=UInt32),
    (n=32, p=[0,19,0],  num=10,  T=UInt32),
    (n=32, p=[0,14,5],  num=10,  T=UInt32),
    (n=32, p=[0,0,19],  num=10,  T=UInt32),
  ]
  println("n   | (p_1, p_2, ...)     | Trials | Rank   | Time      | Bytes allocated | T")
  println("----|---------------------|--------|--------|-----------|-----------------|---------")
  
  for test in tests
    p = copy(test.p)
    (time, bytes, _, rank) = experiment(kmc, test)
    println("$(rpad(test.n, 4, " "))| $(rpad(p, 20, " "))| $(rpad(test.num, 7, " "))| $(rpad(rank, 7, " "))| $(rpad(time, 10, " "))| $(rpad(Base.format_bytes(bytes), 16))| $(test.T)")
  end
end





function benchmark_erect()
  tests = [
    (n=8,  p=[0,4,0],   T=UInt16),
    (n=10, p=[0,6,0],   T=UInt16),
    (n=10, p=[0,4,2],   T=UInt16),
    (n=10, p=[0,0,6],   T=UInt16),
    (n=10, p=[0,1,1,1], T=UInt16),
    (n=13, p=[0,8,0],   T=UInt16),
    (n=13, p=[0,5,3],   T=UInt16),
    (n=13, p=[0,0,8],   T=UInt16),
    (n=16, p=[0,10,0],  T=UInt16),
    (n=16, p=[0,6,4],   T=UInt16),
    (n=16, p=[0,0,10],  T=UInt16),
    (n=24, p=[0,14,0],  T=UInt16),
    (n=24, p=[0,10,4],  T=UInt16),
    (n=24, p=[0,0,14],  T=UInt16),
    (n=32, p=[0,19,0],  T=UInt32),
    (n=32, p=[0,14,5],  T=UInt32),
    (n=32, p=[0,0,19],  T=UInt32),
  ]
    
  println("n   | (p_1, p_2, ...)     | Trials | Time (s)  | Bytes allocated | GC time")
  println("----|---------------------|--------|-----------|-----------------|---------")

  for test in tests 
    p = copy(test.p)
    expr = :( @benchmark my_random_erection($(test.n), $p, $(test.T)) seconds=600 )
    b = eval(expr)

    time = median(b).time
    if time > 1e9
      time = "$(round(time / 1e9, digits=4))s"
    elseif time > 1e6
      time = "$(round(time / 1e6, digits=4))ms"
    elseif time > 1e3
      time = "$(round(time / 1e3, digits=4))µs"
    else
      time = "$(round(time, digits=4))ns"
    end

    bytes = Base.format_bytes(median(b).memory)
    gc = round(median(b).gctime / median(b).time, digits=2)
    trials = length(b.times)


    println("$(rpad(test.n, 4, " "))| $(rpad(p, 20, " "))| $(rpad(trials, 7, " "))| $(rpad(time, 10, " "))| $(rpad(bytes, 16))| $gc%")
  end
end




function generate_benchmark_table_kmc(kmc)
  tests = [
    (n=10, p=[0,6,0], t=100, T=UInt16)
    (n=10, p=[0,0,6], t=100, T=UInt16)
    (n=10, p=[0,1,1,1], t=100, T=UInt16)

    (n=16, p=[0,8,0], t=10, T=UInt16)
    (n=16, p=[0,7,2], t=10, T=UInt16)
    (n=16, p=[0,6,2,2], t=10, T=UInt16)

    (n=32, p=[0,19,0], t=10, T=UInt32)
    (n=32, p=[0,15,4,2], t=10, T=UInt32)
    (n=32, p=[0,0,15,7], t=10, T=UInt32)
    (n=32, p=[0,2,2,2,2,2], t=10, T=UInt32)

    (n=64, p=[0,38,0], t=1, T=UInt64)
    (n=64, p=[0,28,8], t=1, T=UInt64)
    (n=64, p=[0,0,24,8,8], t=1, T=UInt64)

    (n=128, p=[0,77,0], t=1, T=UInt128)
    (n=128, p=[0,56,16], t=1, T=UInt128)
    (n=128, p=[0,0,48,16,16], t=1, T=UInt128)
  ]

  println("n   | (p_1, p_2, ...)     | Trials | Time      | GC time   | Bytes allocated | T")
  println("----|---------------------|--------|-----------|-----------|-----------------|---------")

  for test in tests
    (time, bytes, gctime) = experiment(kmc, test.t, test.n, test.p, test.T)
    println("$(rpad(test.n, 4, " "))| $(rpad(test.p, 20, " "))| $(rpad(test.t, 7, " "))| $(rpad(time, 10, " "))| $(rpad(gctime, 10, " "))| $(rpad(Base.format_bytes(bytes), 16))| $(test.T)")
  end

end

  function generate_benchmark_table_erect(erect)
    tests = [
    (n=10, p=[0,6,0], t=100, T=UInt16)
    (n=10, p=[0,5,1], t=100, T=UInt16)
    (n=10, p=[0,5,2], t=100, T=UInt16)
    (n=10, p=[0,6,1], t=100, T=UInt16)
    (n=10, p=[0,4,2], t=100, T=UInt16)
    (n=10, p=[0,3,3], t=100, T=UInt16)
    (n=10, p=[0,0,6], t=100, T=UInt16)
    (n=10, p=[0,1,1,1], t=100, T=UInt16)
    (n=13, p=[0,6,0], t=100, T=UInt16)
    (n=13, p=[0,6,2], t=100, T=UInt16)
    (n=16, p=[0,6,0], t=100, T=UInt16)
    (n=16, p=[0,6,1], t=100, T=UInt16)
    (n=16, p=[0,0,6], t=100, T=UInt16)
    (n=20, p=[0,6,0], t=10, T=UInt32)
    (n=20, p=[0,6,2], t=10, T=UInt32)
    (n=32, p=[0,6,2,1], t=10, T=UInt32)
    (n=64, p=[0,6,4,4,2,1], t=1, T=UInt64)
    (n=128, p=[0,6,6,4,4,2,1], t=1, T=UInt128)
    (n=129, p=[0,6,6,4,4,2,1], t=1, T=BigInt)
  ]

  println("n   | (p_1, p_2, ...)      | Trials | Rank | Bases | Circuits | Time      | GC time   | Bytes allocated")
  println("----|----------------------|--------|------|-------|----------|-----------|-----------|-----------------")

  for test in tests
    matroids = Dict()
    for _ in 1:test.t
      res = @timed m = erect(test.n, test.p, test.T)
      rank = length(m.F)-1
      if !haskey(matroids, rank) 
        matroids[rank] = Dict("trials"=>0, "bases"=>[], "circuits"=>[], "time"=>[], "gctime"=>[], "bytes"=>[])
      end

      matroids[rank]["trials"] += 1
      push!(matroids[rank]["bases"], length(m.I[rank+1]))
      push!(matroids[rank]["circuits"], length(m.C))
      push!(matroids[rank]["time"], res.time)
      push!(matroids[rank]["gctime"], res.gctime)
      push!(matroids[rank]["bytes"], res.bytes)
    end

    for rank in sort(collect(keys(matroids)))
      matroids[rank]["bases"] = round(avg(matroids[rank]["bases"]), digits=1)
      matroids[rank]["circuits"] = round(avg(matroids[rank]["circuits"]), digits=1)
      matroids[rank]["time"] = round(avg(matroids[rank]["time"]), digits=7)
      matroids[rank]["gctime"] = round(avg(matroids[rank]["gctime"]), digits=7)
      matroids[rank]["bytes"] = round(avg(matroids[rank]["bytes"]), digits=7)

      print("$(rpad(test.n, 4))| $(rpad(test.p, 21))| ")
      print("$(rpad(matroids[rank]["trials"], 7))| ")
      print("$(rpad(rank, 5))| ")
      print("$(rpad(matroids[rank]["bases"], 6))| ")
      print("$(rpad(matroids[rank]["circuits"], 9))| ")
      print("$(rpad(matroids[rank]["time"], 10, " "))| ")
      print("$(rpad(matroids[rank]["gctime"], 10, " "))| ")
      print("$(rpad(Base.format_bytes(matroids[rank]["bytes"]), 16))\n")
    end
  end
end