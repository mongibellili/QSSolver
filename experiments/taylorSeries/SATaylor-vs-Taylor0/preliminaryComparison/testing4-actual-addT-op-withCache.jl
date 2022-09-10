import Base: getindex, setindex! 
import Base.:+
#########################################################Normal Taylor################################################
struct Taylor0
    coeffs :: Array{Float64,1}
    order :: Int
    function Taylor0(coeffs::Array{Float64,1}, order::Int) 
        #resize_coeffs1!(coeffs, order)
        return new(coeffs, order)
    end
end
#= Taylor0(x::Taylor0)   = x
Taylor0(coeffs::Array{Float64,1})  = Taylor0(coeffs, length(coeffs)-1)
function Taylor0(x::Float64, order::Int) 
    v = fill(zero(x), order+1)
    v[1] = x
    return Taylor0(v, order)
end
Taylor0(::Float64, order::Int)  = Taylor0( [0.0, 1.0], order)
Taylor0(order::Int) = Taylor0(Float64, order) =#
getindex(a::Taylor0, n::Int64) = a.coeffs[n+1]
setindex!(a::Taylor0, x::Float64, n::Int)  = a.coeffs[n+1] = x

function (+)(a::Taylor0, b::T) where {T<:Number}
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = (+)(a[0], b)
    return Taylor0(coeffs, a.order)
end

################################################################SA taylor##########################################
using StaticArrays
struct TaylorSA{V} #<: AbstractSeries
    coeffs :: MVector{V,Float64}
    order :: Int
    function TaylorSA(coeffs::MVector{V,Float64}, order::Int) where {V}
        #resize_coeffs1!(coeffs, order)
        return new{V}(coeffs, order)
    end
end
#= TaylorSA(x::TaylorSA)  = x
TaylorSA(coeffs::MVector{V,Float64}) where {V} = TaylorSA(coeffs, length(coeffs)-1)
#= function TaylorSA(x::Float64, order::Int) 
    #MVector{order+1,Int}(vect)
    v=@MVector zeros(order+1)
    #v = fill(zero(x), order+1)
    v[1] = x
    return TaylorSA(v, order)
end =#
TaylorSA(order::Int) = TaylorSA(1.0, order) =#
getindex(a::TaylorSA, n::Int64) = a.coeffs[n+1]
setindex!(a::TaylorSA, x::Float64, n::Int)  = a.coeffs[n+1] = x
#= tSA2= TaylorSA(2.3,5)
tSA3= TaylorSA(5)  =#
#= function (+)(a::TaylorSA, b::Float64) 
    #coeffs = copy(a.coeffs) # removing this will change also a !!!?? a.coeffs passed by reference. copy is needed, so there must be a cache.
    @inbounds a.coeffs[1] = (+)(a[0], b)
    return TaylorSA(a.coeffs, a.order)
end =#
function (+)(a::TaylorSA, b::Float64) 
    coeffs = copy(a.coeffs) # removing this will change also a !!!?? a.coeffs passed by reference. copy is needed, so there must be a cache.
    @inbounds coeffs[1] = (+)(a[0], b)
    return TaylorSA(coeffs, a.order)
end
############################################################testing#########################################"
using BenchmarkTools
using InteractiveUtils
vect=(1.2,2.2,3.0)
coeffs=[vect...]
coeffs2=[1.0,2.0,3.0]
function addT(a::Taylor0, b::Taylor0,c::T,cache::Taylor0) where {T<:Number}
    #cache.coeffs.=a.coeffs  #not needed for in middle ops
  #=   cache.coeffs[1]=a.coeffs[1]
    cache.coeffs[2]=a.coeffs[2]
    cache.coeffs[3]=a.coeffs[3]
    cache.coeffs[1]=cache.coeffs[1]+c  =#
    #@__dot__ cache.coeffs = (+)(cache.coeffs, b.coeffs)
    cache.coeffs[1] = a.coeffs[1]+b.coeffs[1]+c
    cache.coeffs[2] = a.coeffs[2]+b.coeffs[2]
    cache.coeffs[3] = a.coeffs[3]+b.coeffs[3]
      return cache
end  
function subT(a::Taylor0, b::Taylor0,cache::Taylor0) where {T<:Number}
   # @__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
    cache.coeffs[1] = a.coeffs[1]-b.coeffs[1]
    cache.coeffs[2] = a.coeffs[2]-b.coeffs[2]
    cache.coeffs[3] = a.coeffs[3]-b.coeffs[3]
      return cache
end  
t1=Taylor0(coeffs,2)
t2=Taylor0(coeffs2,2)
cache=Taylor0([0.0,0.0,0.0],2)
#= @show addT(t1,t2,1.0,cache)
@btime addT(t1,t2,1.0,cache) =#

############################################################testing#########################################"
staticCoeffs=MVector{3,Float64}(1.2,2.2,3.0)
staticCoeffs2=MVector{3,Float64}(1.0,2.0,3.0)
cachestaticCoeffs2=MVector{3,Float64}(0.0,0.0,0.0)
function addT2(a::TaylorSA, b::TaylorSA,c::T,cache::TaylorSA) where {T<:Number}
    #cache.coeffs.=a.coeffs  #not needed for in middle ops
    cache.coeffs[1]=a.coeffs[1]+b.coeffs[1]+c
    cache.coeffs[2]=a.coeffs[2]+b.coeffs[2]
    cache.coeffs[3]=a.coeffs[3]+b.coeffs[3]
    #cache[0]=cache[0]+c 
   # @__dot__ cache.coeffs = (+)(cache.coeffs, b.coeffs)
      return cache
end  
tSA1= TaylorSA(staticCoeffs,2)
tSA2=TaylorSA(staticCoeffs2,2)
staticcache=TaylorSA(cachestaticCoeffs2,2)


#= @show addT(t1,t2,1.0,cache)
@show addT2(tSA1,tSA2,1.0,staticcache) =#
@btime addT(t1,t2,1.0,cache)#19.822 ns (0 allocations: 0 bytes)
@btime addT2(tSA1,tSA2,1.0,staticcache)#21.995 ns (0 allocations: 0 bytes)
#display(@code_typed subT(t1,t2,cache))
#display(@code_typed addT2(tSA1,tSA2,1.0,staticcache))