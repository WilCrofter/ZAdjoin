import Base: +, *, one, string, show

## Types

""" Modulus

   Assuming α^n = p[1] + p[2]α + ... + p[n]α^(n-1), create expressions for α^k, k=n, ... ,2n-2
"""
immutable Modulus{I<:Integer}
  modulus::Array{Array{I,1},1}
  function (::Type{Modulus}){I}(p::Array{I,1})
    n = length(p)
    n > 1 || error("n==1, a trivial case, is not supported")
    modulus = Array(Array{I,1},n-1)
    modulus[1]=p
    if n>1
      for i in 2:(n-1)
        modulus[i] = vcat([0],modulus[i-1][1:(n-1)])+modulus[i-1][1]*p
      end
    end
    new{I}(modulus)
  end
end

immutable Element{I<:Integer}
  p::Array{I,1}
  modulus::Modulus{I}
  
  function (::Type{Element}){I}(p::Array{I,1}, modulus::Modulus{I})
    length(p)==length(modulus.modulus[1]) || error("unequal degrees ", length(p), " and ", length(modulus.modulus[1]),".")
    new{I}(p,modulus)
  end
end

## Arithmetic

function +{I<:Integer}(a::Element{I},b::Element{I})
  a.modulus == b.modulus || error("terms are not from the same ring")
  Element(a.p+b.p,a.modulus)
end

function *{I<:Integer}(a::Element{I},b::Element{I})
  a.modulus == b.modulus || error("factors are not from the same ring")
  n=length(a.p)
  ans = zeros(I,2*n-1)
  # convolve (noting that arrays are 1-origin)
  for k in 2:(2*n)
     for i in max(k-n,1):min(k-1,n)
       ans[k-1]+=a.p[i]*b.p[k-i]
     end
  end
  for i in 1:(n-1)
    ans[1:n] += ans[n+i]*a.modulus.modulus[i]
  end
  Element(ans[1:n], a.modulus)
end

function one{I<:Integer}(z::Element{I})
  p=zeros(I,length(z.p))
  p[1]=1
  Element(p,z.modulus)
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

