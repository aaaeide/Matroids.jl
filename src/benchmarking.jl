function avg(arr)
  return sum(arr) / length(arr)
end

function benchmark(kmc, num, n, p, T)
  times = []
  bytes = []
  gctimes = []
  for _ in 1:num
    res = @timed kmc(n, p, T)
    push!(times, res.time)
    push!(bytes, res.bytes)
    push!(gctimes, res.gctime)
  end

  t = round(avg(times), digits=7)
  b = floor(Int, avg(bytes))
  g = round(avg(gctimes), digits=7)

  return (t,b,g)
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
    (time, bytes, gctime) = benchmark(kmc, test.t, test.n, test.p, test.T)
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