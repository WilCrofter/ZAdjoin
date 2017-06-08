## Generator

function generator{I<:Integer}(modulus::AbstractModulus{I})
  p = zeros(I,length(modulus.polynomials[1]))
  p[2]=I(1)
  Element(p, modulus)
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

function det{I<:Integer}(z::Element{I})
  I(round(det(Matrix(z))))
end

## Conveniences

function +{I<:Integer, J<:Integer}(a::Element{I}, b::J)
  p=deepcopy(a.coefficients)
  p[1] += I(b)
  Element(p,a.modulus)
end

+{I<:Integer, J<:Integer}(b::J, a::Element{I}) = a+I(b)
*{I<:Integer, J<:Integer}(a::Element{I}, b::J) = Element(I(b)*deepcopy(a.coefficients), a.modulus)
*{I<:Integer, J<:Integer}(b::J, a::Element{I}) = a*I(b)

## Equality

function =={I<:Integer}(a::Element{I}, b::Element{I})
  a.coeffficients == b.coefficients && a.modulus == b.modulus
end
