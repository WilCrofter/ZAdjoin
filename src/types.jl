""" Modulus

  Assuming α^n = p[1] + p[2]α + ... + p[n]α^(n-1), create expressions for α^k, k=n, ... ,2n-2 and powers 0,...,n-1of p's companion matrix.
  """
abstract AbstractModulus{I<:Integer}

immutable Modulus{I<:Integer} <: AbstractModulus{I}
  polynomials::Array{Array{I,1},1}
  matrices::Array{Array{I,2},1}
  
  function (::Type{Modulus}){I}(p::Array{I,1})
    n = length(p)
    n > 1 || error("n==1, a trivial case, is not supported")
    polynomials,matrices = modulusfields(p)
    new{I}(polynomials,matrices)
  end
  
end

immutable KummerModulus{I<:Integer} <: AbstractModulus{I}
  polynomials::Array{Array{I,1},1}
  matrices::Array{Array{I,2},1}
  
  function (::Type{KummerModulus}){I}(n::I)
    k = n%2==0 ? Int(n/2) : n
    k > 1 || error("n==$n, a trivial case, is not supported")
    i = n%2==0 ? I(-1) : I(1)
    polynomials,matrices = modulusfields(vcat(i,zeros(I,k-1)))
    new{I}(polynomials,matrices)
  end
  
end

function modulusfields{I}(p::Array{I,1})
  n=length(p)
  polynomials = Array(Array{I,1},n-1)
  polynomials[1]=p
  A = zeros(I,(n,n))
  A[:,n]=p
  if n>1
    for i in 2:(n-1)
      polynomials[i] = vcat([0],polynomials[i-1][1:(n-1)])+(polynomials[i-1][n])*p
      A[i,i-1]=1
    end
    A[n,n-1]=1
  end
  matrices = Array(Array{I,2},n)
  for i in 1:n
    matrices[i]=A^(i-1)
  end
  return polynomials, matrices
end

immutable Element{I<:Integer}
  coefficients::Array{I,1}
  modulus::AbstractModulus{I}
  
  function (::Type{Element}){I}(coef::Array{I,1}, modulus::AbstractModulus{I})
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

function -{I<:Integer}(a::Element{I},b::Element{I})
  a.modulus == b.modulus || error("terms are not from the same ring")
  Element(a.coefficients-b.coefficients,a.modulus)
end

function +{I<:Integer}(a::Element{I})
  Element(+(a.coefficients), a.modulus)
end

function -{I<:Integer}(a::Element{I})
  Element(-(a.coefficients), a.modulus)
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

function string{I<:Integer}(m::AbstractModulus{I})
  "α^$(length(m.polynomials[1])) = $(helper(m.polynomials[1]))"
end

function show{I<:Integer}(io::IO, z::Element{I})
  print(io, "$(typeof(z))\n")
  print(io, string(z))
end

function show{I<:Integer}(io::IO, m::AbstractModulus{I})
  print(io, "$(typeof(m))\n")
  print(io, string(m))
end
