import Base: +, *, one, string, show, Matrix

## Types

""" Modulus

   Assuming α^n = p[1] + p[2]α + ... + p[n]α^(n-1), create expressions for α^k, k=n, ... ,2n-2
"""
immutable Modulus{I<:Integer}
  modulus::Array{Array{I,1},1}
  matrix::Array{I,2}
  
  function (::Type{Modulus}){I}(p::Array{I,1})
    n = length(p)
    n > 1 || error("n==1, a trivial case, is not supported")
    modulus = Array(Array{I,1},n-1)
    matrix = zeros(I,(n,n))
    modulus[1]=p
    matrix[:,n]=p
    if n>1
      for i in 2:(n-1)
        modulus[i] = vcat([0],modulus[i-1][1:(n-1)])+modulus[i-1][1]*p
        matrix[i,i-1]=1
      end
      matrix[n,n-1]=1
    end
    new{I}(modulus,matrix)
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

## Conversion

function Matrix{I}(z::Element{I})
  p = z.p
  n = length(p)
  A = z.modulus.matrix
  B = eye(I,n)
  ans = p[1]*B
  for i=2:n
    B*=A
    ans += p[i]*B
  end
  ans
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

function helper{I<:Integer}(p::Array{I,1})
  tmp=""
  for i in eachindex(p)
    if p[i] != 0
      if i==1
        tmp=tmp*"$(p[1])"
      elseif i==2
        tmp=tmp*"+$(p[i])α"
      else
        tmp=tmp*"+$(p[i])α^$(i-1)"
      end
    end
  end
  return length(tmp) == 0 ? "0" : tmp
end

function string{I<:Integer}(z::Element{I})
  helper(z.p)
end

function string{I<:Integer}(m::Modulus{I})
  "α^$(length(m.modulus[1])) = $(helper(m.modulus[1]))"
end

function show{I<:Integer}(io::IO, z::Element{I})
  print(io, "Element{$I}\n")
  print(io, string(z))
end

function show{I<:Integer}(io::IO, m::Modulus{I})
  print(io, "Modulus{$I}\n")
  print(io, string(m))
end

