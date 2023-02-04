
"""
Second implementation of Knuth's matroid construction algorithm, this
time faithful to Knuth's ERECTION.W.
"""
function erection(n, enlargements)
  avail = nothing # Beginning of list of available space
  nh = nothing # Head of circular list being formed for rank |r+1|

  rank = initialize_rank(n)
  L, S, unused = initialize_list_memory(n)

  rank[1] = 0; r = 0
  while rank[1<<n - 1] > r
    # Create empty list
    nh = avail
    if nh != nothing
      avail = L[nh]
    else
      nh = unused + 1
    end
    L[nh] = nh


  end
end

function initialize_rank(n)
  k = 1
  rank = zeros(1<<n)
  rank[1] = 100

  while k < 1<<n
    for i in 0:k-1
      rank[k+i+1] = rank[i+1]+1
    end
    k = k+k
  end

  return rank
end

function initialize_list_memory(n)
  lmax = 25742 # 2((16C8)+1), a safe upper bound on list size
  L = zeros(lmax)
  S = zeros(lmax)
  L[1]=2; L[2]=1; S[2]=0; h=1

  return L, S, 3
end
