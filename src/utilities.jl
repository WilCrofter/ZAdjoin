function issquarefree{I<:Integer}(d::I)
  # factor(d) is a dictionary whose keys are primes
  # and whose values are corresponding prime powers.
  # An integer is square free iff none of its prime
  # factors have power 2 or greater. This is equivalent
  # to its Mobius function being nonzero.
  mobius(d) != 0
end

"""

[Link](https://proofwiki.org/wiki/Definition:M%C3%B6bius_Function)
"""
function mobius{I<:Integer}(n::I)
  if n==1 return I(1) end
  powers = collect(values(factor(n)))
  return all(powers.== 1) ? I(-1)^length(powers) : I(0)
end


#  Division by monic polynomials

function onestep!{I<:Integer}(num::Array{I,1}, den::Array{I,1}, quo::Array{I,1})
  n = findlast(num .!= 0)
  d = findlast(den .!= 0)
  den[d] == 1 || error("denominator is not monic")
  k = n-d
  k >= 0 || error("numerator not divisible by denominator")
  quo[1+k] += num[n]
  for i=1:d
    num[i+k] -= num[n]*den[i]
  end
end

function divide_by_monic{I<:Integer}(num::Array{I,1}, den::Array{I,1}) 
  cnum = deepcopy(num)
  quo = zeros(I,length(num))
  while any(cnum .!= 0)
    onestep!(cnum, den, quo)
  end
  quo
end

function cyclotomic{I<:Integer}(n::I)
  n >= 1 || error("There is no cyclotomic of degree < 1")
  if n==1
    return [I(-1), I(1)]
  end
  num = zeros(I,n+1)
  num[n+1]=I(1)
  num[1] = I(-1)
  for d in 1:(n-1)
    if n%d == 0
      num = divide_by_monic(num, cyclotomic(d))
    end
  end
  d = findlast(num .!= 0)
  num[1:d]
end
