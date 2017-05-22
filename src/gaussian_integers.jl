#= Use Gaussian integers as an example to work out type implemenations.
=#

import Base.size, Base.getindex, Base.setindex!, Base.ndims, Base.print_matrix
import Base.+

immutable GInt <: AbstractArray{Integer,2}
  data::Array{Integer,2}
  function GInt{I<:Integer}(data::Array{I,2})
    new(data)
  end
end

function +(x::GInt,y::GInt)
  tmp = x.data+y.data
  GInt(tmp[1,1],tmp[1,2])
end

function print_matrix(ioc::IOContext{Base.Terminals.TTYTerminal}, g::RingSandbox.GInt, s1::String, s2::String, s3::String)
  show(ioc,"$(g.data[1,1])I + $(g.data[1,2])α")
end

function ndims(g::GInt)
  ndims(g.data)
end

function size(g::GInt)
  size(g.data)
end

function getindex(g::GInt, i::Int, j::Int)
  0<i<3 && 0<j<3 || error("illegal index")
  g.data[i,j]
end

function getindex(g::GInt, I...)
  getindex(g.data,I)
end

function setindex!(g::GInt, v, i::Int)
  error("GInts are immutable")
end

function setindex!(g::GInt, v, I...)
  error("GInts are immutable")
end

