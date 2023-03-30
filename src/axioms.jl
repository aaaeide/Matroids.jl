include("kmc.jl")

"""
The ground set is closed. E ∈ F.
"""
function c1(m::KnuthMatroid)
  @assert last(m.F) == Set(big"2"^m.n-1)
end

"""
The intersection of two closed sets is a closed set. 
If A, B ∈ F, then A ∩ B ∈ F.
"""
function c2(m::KnuthMatroid)
  F = reduce(∪, m.F)
  for A ∈ F for B ∈ F
    @assert A & B ∈ F "$A ∩ $B = ∉ F"
  end end
end

"""
If A ∈ F and a, b ∈ E - A, then b is a member of all sets containing A ∪ {a} if and only if a is a member of all sets containing A ∪ {b}.
"""
function c3(m::KnuthMatroid)
  E = m.Type(big"2"^m.n-1)
  F = reduce(∪, m.F)
  for A ∈ setdiff(F, 0)
    t1 = E-A
    while t1 > 0
      a = t1&-t1
      t2 = t1&~a
      while t2 > 0
        b = t2&-t2
        ā = reduce(&, [B for B in F if (A|a)&B == A|a])
        b̄ = reduce(&, [B for B in F if (A|b)&B == A|b])
        @assert (b&ā==b) == (a&b̄==a)
        t2 &= ~b
      end
      t1 &= ~a
    end
  end
end