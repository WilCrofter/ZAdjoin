"""
    ZAdjoin

  Let α be an element satisfying a monic polynomial of degree n > 1 over the integers, Z. That is,

  ```
  α^n = p[1] + p[2]α + ... + p[n]α^(n-1),  (I)
  ```

  where the p[i]ϵZ. Let Z[α] be the ring formed by adjoining α to the Z. This module provides code for basic representation and manipulation of such rings.

  Elements of Z[α] can be represented as polynomials, q[1] + q[2]α + ... + q[n]α^(n-1), of degree n-1 in α, since higher powers of α can be reduced modulo the defining relationship, α^n=p(α). Thus, under addition, Z[α] is a module of dimension n over Z. Since, under multiplication, elements of Z[α] induce linear transformations of this module, elements may also be represented as matrices over Z. Specifically, α may be represented as the companion matrix of the defining polynomial, α^n-p(α), powers of α up to n-1 as corresponding powers of this companion matrix, and general elements as integer combinations of these powers.

  A defining relationship, essentially n integer coefficients p[1],...,p[n] as in equation (I) above, is represented by one of the types Modulus, KummerModulus, or QIModulus (QI for Quadratic Integer,) all concrete subtypes of AbstractModulus. Ring elements are represented as integer vectors and a modulus (instance of a concrete subtype of AbstractModulus) combined in the Element type. Addition, subtraction, multiplication, and exponentiation by integer powers of Elements are supported. Conversion of Elements to matrices is also supported.

  # Examples

  ```
  # A Kummer Ring of degree 4 is the Gaussian (complex) Integers
  julia> g = KummerModulus(4) 
  ZAdjoin.KummerModulus{Int64}
  α^2 = -1

  # Function generator(modulus) will return the adjoined element.
  julia> α = generator(g)

  # Other elements can be formed from the adjoined element.
  julia> β = 3+4α # Another element
  ZAdjoin.Element{Int64}
  3+4α

  julia> α + β
  ZAdjoin.Element{Int64}
  3+5α

  julia> α*β
  ZAdjoin.Element{Int64}
  -4+3α

  # Compare with complex multiplication
  julia> im*(3+4im)
  -4 + 3im

  julia> α^2
  ZAdjoin.Element{Int64}
  -1

  # Matrix representation
  julia> Matrix(β)
  2×2 Array{Int64,2}:
   3  -4
   4   3

  # Representations are arithmetcally compatible (isomorphic)
  julia> Matrix(α+β^5)
  2×2 Array{Int64,2}:
    -237  3115
   -3115  -237

  julia> Matrix(α)+Matrix(β)^5
  2×2 Array{Int64,2}:
    -237  3115
   -3115  -237

  # The determinant of an element is the determinant of its Matrix
  # and acts as a multiplicative homomorphism into the integers.
  julia> det(β)
  25
  julia> det(1+α)
  2
  julia> det((1+α)*β)
  50

  ```


"""
module ZAdjoin

import Base: +, -, *, one, det, ==, string, show, Matrix
import Primes
export AbstractModulus, Modulus, KummerModulus, Element, generator

include("types.jl")
include("operations.jl")
include("utilities.jl")


end
