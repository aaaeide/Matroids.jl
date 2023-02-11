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

  println("time\t$t\ngctime\t$g\nbytes\t$b")
end