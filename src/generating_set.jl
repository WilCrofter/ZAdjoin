

immutable Generator{I<:Integer}
  data::Array{Array{I,2}}
  
  function (::Type{Generator}){I}(p::Array{I,1})
    N=length(p)
    X = zeros(I,(N,N))
    d = Array(Array{I,2},N)
    d[1]=eye(I,N)
    if N > 1
      X[:,1]=-p
      for i in 2:N
        X[i-1,i]=1
      end
      for i in 2:N
        d[i]=d[i-1]*X
      end
    end
    new{I}(d)
  end
  
end

