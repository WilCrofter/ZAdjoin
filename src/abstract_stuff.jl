import Base.size, Base.getindex, Base.setindex!, Base.ndims, Base.print_matrix
import Base.+

immutable SquaresVector <: AbstractArray{Int, 1}
  count::Int
end

Base.size(S::SquaresVector) = (S.count,)
Base.linearindexing{T<:SquaresVector}(::Type{T}) = Base.LinearFast()
Base.getindex(S::SquaresVector, i::Int) = i*i;

immutable GI <: AbstractArray{Int, 2}
  data::Array{Int,2}
  GI(a::Array{Int,2}) = size(a)==(2,2) ? new(a) : error("wrong size")
end


Base.size(g::GI)=size(g.data)
Base.linearindexing{T<:GI}(::Type{T})=Base.LinearFast()
Base.getindex(g::GI,i::Int,j::Int) = g.data[i,j]


immutable SparseArray{T,N} <: AbstractArray{T,N}
  data::Dict{NTuple{N,Int}, T}
  dims::NTuple{N,Int}
end

SparseArray{T}(::Type{T}, dims::Int...) = SparseArray(T, dims)
SparseArray{T,N}(::Type{T}, dims::NTuple{N,Int}) = SparseArray{T,N}(Dict{NTuple{N,Int}, T}(), dims)

Base.size(A::SparseArray) = A.dims
Base.similar{T}(A::SparseArray, ::Type{T}, dims::Dims) = SparseArray(T, dims)
# Define scalar indexing and indexed assignment
Base.getindex{T,N}(A::SparseArray{T,N}, I::Vararg{Int,N})     = get(A.data, I, zero(T))
Base.setindex!{T,N}(A::SparseArray{T,N}, v, I::Vararg{Int,N}) = (A.data[I] = v)
