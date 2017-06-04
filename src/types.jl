""" Modulus

  Assuming α^n = p[1] + p[2]α + ... + p[n]α^(n-1), create expressions for α^k, k=n, ... ,2n-2 and powers 0,...,n-1of p's companion matrix.
  """
immutable Modulus{I<:Integer}
  polynomials::Array{Array{I,1},1}
  matrices::Array{Array{I,2},1}
  
  function (::Type{Modulus}){I}(p::Array{I,1})
    n = length(p)
    n > 1 || error("n==1, a trivial case, is not supported")
    polynomials = Array(Array{I,1},n-1)
    polynomials[1]=p
    A = zeros(I,(n,n))
    A[:,n]=p
    if n>1
      for i in 2:(n-1)
        polynomials[i] = vcat([0],polynomials[i-1][1:(n-1)])+polynomials[i-1][1]*p
        A[i,i-1]=1
      end
      A[n,n-1]=1
    end
    matrices = Array(Array{I,2},n)
    for i in 1:n
      matrices[i]=A^(i-1)
    end
    new{I}(polynomials,matrices)
  end
  
end

immutable Element{I<:Integer}
  coefficients::Array{I,1}
  modulus::Modulus{I}
  
  function (::Type{Element}){I}(coef::Array{I,1}, modulus::Modulus{I})
    length(coef)==length(modulus.polynomials[1]) || error("unequal degrees ", length(coef), " and ", length(modulus.polynomials[1]),".")
    new{I}(coef,modulus)
  end
end

## Conversion

function Matrix{I}(z::Element{I})
  p = z.coefficients
  n = length(p)
  matrices = z.modulus.matrices
  ans = p[1]*matrices[1]
  for i in 2:n
    ans += p[i]*matrices[i]
  end
  ans
end

## Arithmetic

function +{I<:Integer}(a::Element{I},b::Element{I})
  a.modulus == b.modulus || error("terms are not from the same ring")
  Element(a.coefficients+b.coefficients,a.modulus)
end

function *{I<:Integer}(a::Element{I},b::Element{I})
  a.modulus == b.modulus || error("factors are not from the same ring")
  n=length(a.coefficients)
  ans = zeros(I,2*n-1)
  # convolve (noting that arrays are 1-origin)
  for k in 2:(2*n)
     for i in max(k-n,1):min(k-1,n)
       ans[k-1]+=a.coefficients[i]*b.coefficients[k-i]
     end
  end
  for i in 1:(n-1)
    ans[1:n] += ans[n+i]*a.modulus.polynomials[i]
  end
  Element(ans[1:n], a.modulus)
end

function one{I<:Integer}(z::Element{I})
  p=zeros(I,length(z.coefficients))
  p[1]=1
  Element(p,z.modulus)
end

## Equality

function =={I<:Integer}(a::Element{I}, b::Element{I})
  a.coeffficients == b.coefficients && a.modulus == b.modulus
end

## Display

function helper{I<:Integer}(p::Array{I,1})
  tmp=""
  for i in eachindex(p)
    if p[i] != 0
      sgn = p[i] < 0 ? "-" : "+"
      if i==1
        tmp=tmp*"$(p[1])"
      elseif i==2
        tmp=tmp*"$(sgn)$(abs(p[i]))α"
      else
        tmp=tmp*"$(sgn)$(abs(p[i]))α^$(i-1)"
      end
    end
  end
  return length(tmp) == 0 ? "0" : tmp
end

function string{I<:Integer}(z::Element{I})
  helper(z.coefficients)
end

function string{I<:Integer}(m::Modulus{I})
  "α^$(length(m.polynomials[1])) = $(helper(m.polynomials[1]))"
end

function show{I<:Integer}(io::IO, z::Element{I})
  print(io, "Element{$I}\n")
  print(io, string(z))
end

function show{I<:Integer}(io::IO, m::Modulus{I})
  print(io, "Modulus{$I}\n")
  print(io, string(m))
end
