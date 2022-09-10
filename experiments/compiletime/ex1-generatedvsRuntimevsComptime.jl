using CompTime
using BenchmarkTools
struct SVector{T,n}
    v::Vector{T}
    function SVector(v::Vector{T}) where {T}
      new{T,length(v)}(v)
    end
    function SVector{T,n}(v::Vector{T}) where {T,n}
      #assert(n == length(v))
      new{T,n}(v)
    end
    function SVector{T,n}() where {T,n}
      new{T,n}(Vector{T}(undef,n))
    end
  end
function Base.getindex(sv::SVector{Float64, 3}, i::Int64)
    return getindex(sv.v, i)
end
function Base.setindex!(sv::SVector{Float64, 3},f::Float64, i::Int64)
    return setindex!(sv.v,f, i)
end
 
  

  
  @generated function generatedAdd( v1::SVector{T,n}, v2::SVector{T,n}) where {T,n}
      Expr(:block,
      :(vout = SVector{$T,$n}(Vector{$T}(undef, $n))),  #code copied from github is missing $n; I added it
      begin
        code = Expr(:block)
        for i in 1:n
          push!(code.args, :(vout[$i] = v1[$i] + v2[$i]))
        end
        code
      end,
      :(vout)
    )
  end

  function runtime( v1::SVector{T,n}, v2::SVector{T,n}) where {T,n}
    vout = SVector{T,n}()
    for i in 1:n
      vout[i] = v1[i] + v2[i]
    end
    vout
  end

  @ct_enable function comptimeadd(v1::SVector{T,n}, v2::SVector{T,n}) where {T,n}
    vout = SVector{(@ct T), (@ct n)}()
    @ct_ctrl for i in 1:n
      vout[@ct i] = v1[@ct i] + v2[@ct i]
    end
    vout
  end
  
   

  v1=[2.2,1.1,3.0]
  v2=[0.2,1.0,0.0]
sv1=SVector(v1)
sv2=SVector(v2)

@btime runtime(sv1,sv2)#52.784 ns (2 allocations: 96 bytes)
#@show runtime(sv1,sv2)
#@show generatedAdd(sv1,sv2)
@btime generatedAdd(sv1,sv2)#52.784 ns (2 allocations: 96 bytes)
#@show comptimeadd(sv1,sv2)
  @btime comptimeadd(sv1,sv2)#52.751 ns (2 allocations: 96 bytes)