using Base.Test

try
  using ZAdjoin
catch
  info("Testing as a module. Not installed as a package.")
  include("../src/ZAdjoin.jl")
  using ZAdjoin
end

# ZAdjoin.cyclotomic(n) is recursive, requiring cyclotomic
# polynomials for all proper divisors of n (including 1.)
@test ZAdjoin.cyclotomic(24)==[1,0,0,0,-1,0,0,0,1] # 1-x^4+x^8
@test ZAdjoin.cyclotomic(21)==[1,-1,0,1,-1,0,1,0,-1,1,0,-1,1] # 1-x+x^3-x^4+x^6-x^8+x^9-x^11+x^12



nothing

