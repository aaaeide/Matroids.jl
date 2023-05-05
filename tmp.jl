include("src/kmc.jl")


function find_problem(n, p)
  while true
    cs, m1 = random_erect(n, p)
    m2 = random_knuth_matroid(n, p, UInt16, cs)
    if m1.r != m2.r
      print(cs)
      return
    end
  end
end

# Any[Any[], Any[(0x0120, 0x0010)], Any[(0x0260, 0x0004)], Any[(0x0135, 0x0080)]]

function find_discrepancy(n, p, x)
  _, m1 = random_erect(n, p, UInt16, x)
  
  println("---SKILLETEGN---")
  
  m2 = random_knuth_matroid(n, p, UInt16, x)
end


find_problem(10, (0,1,1,1))

# open("out.log", "w") do out
#   redirect_stdout(out) do
#     find_discrepancy(10, (0,1,1,1), Any[Any[], Any[(0x0120, 0x0010)], Any[(0x0260, 0x0004)], Any[(0x0135, 0x0080)]])
#   end
# end