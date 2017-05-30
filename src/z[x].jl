import Base: +, *, string, show

## Types

""" PolyTable

   Assuming x^n = p[1] + p[2]x + ... + p[n]x^(n-1), create expressions for x^k, k=n, ... ,2n-2
"""
immutable PolyTable{I<:Integer}
  tbl::Array{Array{I,1},1}
  function (::Type{PolyTable}){I}(p::Array{I,1})
    n = length(p)
    n > 1 || error("n==1, a trivial case, is not supported")
    tbl = Array(Array{I,1},n-1)
    tbl[1]=p
    if n>1
      for i in 2:(n-1)
        tbl[i] = vcat([0],tbl[i-1][1:(n-1)])+tbl[i-1][1]*p
      end
    end
    new{I}(tbl)
  end
end

immutable Element{I<:Integer}
  p::Array{I,1}
  pt::PolyTable{I}
  
  function (::Type{Element}){I}(p::Array{I,1}, pt::PolyTable{I})
    length(p)==length(pt.tbl[1]) || error("unequal degrees ", length(p), " and ", length(pt.tbl[1]),".")
    new{I}(p,pt)
  end
end

## Arithmetic

function +{I<:Integer}(a::Element{I},b::Element{I})
  a.pt == b.pt || error("terms are not from the same ring")
  Element(a.p+b.p,a.pt)
end

function *{I<:Integer}(a::Element{I},b::Element{I})
  a.pt == b.pt || error("factors are not from the same ring")
  n=length(a.p)
  ans = zeros(I,n)
  for i in eachindex(a.p)
    for j in i:n
      ans[j] += a.p[i]*b.p[j]
    end
    if i > 1
      for j in (n+2-i):n
        ans += a.p[i]*b.p[j]*a.pt.tbl[i+j-n-1] 
      end
    end
  end
  Element(ans, a.pt)
end

## Display

function string{I<:Integer}(z::Element{I})
  tmp=""
  for i in eachindex(z.p)
    if i==1
      tmp=tmp*"$(z.p[1])"
    elseif i==2
      tmp=tmp*"+$(z.p[i])α"
    else
      tmp=tmp*"+$(z.p[i])α^$(i-1)"
    end
  end
  tmp
end

function show{I<:Integer}(io::IO, z::Element{I})
  print(io, string(z))
end

