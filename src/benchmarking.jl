function avg(arr)
  return sum(arr) / length(arr)
end

function benchmark(kmc, num, n, p)
  times = []
  bytes = []
  gctimes = []
  for _ in 1:num
    res = @timed kmc(n, p)
    push!(times, res.time)
    push!(bytes, res.bytes)
    push!(gctimes, res.gctime)
  end

  t = round(avg(times), digits=7)
  b = floor(Int, avg(bytes))
  g = round(avg(gctimes), digits=7)

  return (t,b,g)
end

function generate_benchmark_table(kmc)
  tests = [
    (n=10, p=[0,6,0], t=100)
    (n=10, p=[0,5,1], t=100)
    (n=10, p=[0,5,2], t=100)
    (n=10, p=[0,6,1], t=100)
    (n=10, p=[0,4,2], t=100)
    (n=10, p=[0,3,3], t=100)
    (n=10, p=[0,0,6], t=100)
    (n=10, p=[0,1,1,1], t=100)
    (n=13, p=[0,6,0], t=10)
    (n=13, p=[0,6,2], t=10)
    (n=16, p=[0,6,0], t=10)
  ]

  println("n  | (p_1, p_2, ...) | Trials | Time      | GC time   | Bytes allocated")
  println("---|-----------------|--------|-----------|-----------|-----------------")

  for test in tests
    (time, bytes, gctime) = benchmark(kmc, test.t, test.n, test.p)
    println("$(test.n) | $(rpad(test.p, 16, " "))| $(rpad(test.t, 7, " "))| $(rpad(time, 10, " "))| $(rpad(gctime, 10, " "))| $(Base.format_bytes(bytes))")
  end
end