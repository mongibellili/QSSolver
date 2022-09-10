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

function (+)(a::Taylor0, b::Float64) 
    coeffs = copy(a.coeffs)
    @inbounds coeffs[1] = (+)(a[0], b)
    return Taylor0(coeffs, a.order)
end

################################################################SA taylor##########################################
using StaticArrays
struct TaylorSA #<: AbstractSeries
    coeffs :: MVector{3,Float64}
    order :: Int
    function TaylorSA(coeffs::MVector{3,Float64}, order::Int) 
        #resize_coeffs1!(coeffs, order)
        return new(coeffs, 2)
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
    return TaylorSA(coeffs, 2)
end
############################################################testing#########################################"
using BenchmarkTools
vect=(1.2,2.2,3.0)
coeffs=[vect...]
coeffs2=[1.0,2.0,3.0]
function addT(a::Taylor0, b::Taylor0,c::Float64,cache::Taylor0) 
    cache.coeffs.=a.coeffs  #not needed for in middle ops
      cache[0]=cache[0]+c 
    @__dot__ cache.coeffs = (+)(cache.coeffs, b.coeffs)
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
function addT2(a::TaylorSA, b::TaylorSA,c::Float64,cache::TaylorSA) 
    cache.coeffs.=a.coeffs  #not needed for in middle ops
      cache[0]=cache[0]+c 
    @__dot__ cache.coeffs = (+)(cache.coeffs, b.coeffs)
      return cache
end  
tSA1= TaylorSA(staticCoeffs,2)
tSA2=TaylorSA(staticCoeffs2,2)
staticcache=TaylorSA(cachestaticCoeffs2,2)

function test(t1::Taylor0, t2::Taylor0,c::Float64,cache::Taylor0) 
    addT(t1,t2,c,cache)
end
function testSA(tSA1::TaylorSA, tSA2::TaylorSA,c::Float64,staticcache::TaylorSA) 
    addT2(tSA1,tSA2,c,staticcache)
end

@btime test(t1,t2,1.0,cache)#19.821 ns (0 allocations: 0 bytes)
@btime testSA(tSA1,tSA2,1.0,staticcache)#30.498 ns (1 allocation: 32 bytes)