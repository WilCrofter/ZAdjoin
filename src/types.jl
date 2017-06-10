"""
    AbstractModulus

  

  Assuming α^n = p[1] + p[2]α + ... + p[n]α^(n-1), create expressions for α^k, k=n, ... ,2n-2 and powers 0,...,n-1of p's companion matrix.
  """
abstract AbstractModulus{I<:Integer}

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

immutable Modulus{I<:Integer} <: AbstractModulus{I}
  polynomials::Array{Array{I,1},1}
  matrices::Array{Array{I,2},1}
  symbol::String
  
  function (::Type{Modulus}){I}(p::Array{I,1}; symbol::String="α")
    n = length(p)
    n > 1 || error("n==1, a trivial case, is not supported")
    polynomials,matrices = modulusfields(p)
    new{I}(polynomials,matrices,symbol)
  end
  
end

immutable KummerModulus{I<:Integer} <: AbstractModulus{I}
  polynomials::Array{Array{I,1},1}
  matrices::Array{Array{I,2},1}
  symbol::String
  
  function (::Type{KummerModulus}){I}(n::I; symbol::String="ζ")
    n > 1 || error("n==$n, a trivial case, is not supported")
    p = cyclotomic(n)
    polynomials,matrices = modulusfields(-p[1:(end-1)])
    new{I}(polynomials,matrices,symbol)
  end
  
end

immutable QIModulus{I<:Integer} <: AbstractModulus{I}
  polynomials::Array{Array{I,1},1}
  matrices::Array{Array{I,2},1}
  d::I
  k::I

  function (::Type{QIModulus}){I}(d::I)
    issquarefree(d) || error("$d is not square free")
    if mod(d,4) == 1
      # adjoin (1+sqrt(d))/2
      k = I((d-1)/4)
      polynomials, matrices = modulusfields([k,d])
    else
      # mod(d,4) is either 2 or 3, adjoin sqrt(d)
      polynomials, matrices = modulusfields([d,I(0)])
      k = I(0)
    end
    new{I}(polynomials, matrices, d, k)
  end

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


## Display

function helper{I<:Integer}(p::Array{I,1}, symbol::String)
  tmp=""
  for i in eachindex(p)
    if p[i] != 0
      sgn = p[i] < 0 ? "-" : "+"
      if i==1
        tmp=tmp*"$(p[1])"
      elseif i==2
        tmp=tmp*"$(sgn)$(abs(p[i]))$(symbol)"
      else
        tmp=tmp*"$(sgn)$(abs(p[i]))$(symbol)^$(i-1)"
      end
    end
  end
  return length(tmp) == 0 ? "0" : tmp
end

function string{I<:Integer}(z::Element{I})
  helper(z.coefficients, z.modulus.symbol)
end

function string{I<:Integer}(m::AbstractModulus{I})
  "$(m.symbol)^$(length(m.polynomials[1])) = $(helper(m.polynomials[1],m.symbol))"
end

function show{I<:Integer}(io::IO, z::Element{I})
  print(io, "$(typeof(z))\n")
  print(io, string(z))
end

function show{I<:Integer}(io::IO, m::AbstractModulus{I})
  print(io, "$(typeof(m))\n")
  print(io, string(m))
end
