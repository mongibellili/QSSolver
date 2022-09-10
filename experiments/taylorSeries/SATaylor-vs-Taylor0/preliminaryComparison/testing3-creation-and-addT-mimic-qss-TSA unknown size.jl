import Base: getindex, setindex! 
import Base.:+
#########################################################Normal Taylor################################################
struct Taylor0
    coeffs :: Array{Float64,1}
    order :: Int
    function Taylor0(coeffs::Array{Float64,1}, order::Int) 
        return new(coeffs, order)
    end
end
getindex(a::Taylor0, n::Int64) = a.coeffs[n+1]
setindex!(a::Taylor0, x::Float64, n::Int)  = a.coeffs[n+1] = x
function addT(a::Taylor0, b::Taylor0,c::Float64,cache::Taylor0) #where {T<:Number}
    cache.coeffs[1] = a.coeffs[1]+b.coeffs[1]+c
    cache.coeffs[2] = a.coeffs[2]+b.coeffs[2]
    cache.coeffs[3] = a.coeffs[3]+b.coeffs[3]
      return cache
end
################################################################SA taylor##########################################
using StaticArrays
#= struct TaylorSA #<: AbstractSeries
    coeffs :: MVector{3,Float64}
    order :: Int
    function TaylorSA(coeffs::MVector{3,Float64}, order::Int)
        return new(coeffs, 2)
    end
end =#
struct TaylorSA{V} #<: AbstractSeries
    coeffs :: MVector{V,Float64}
    order :: Int
    function TaylorSA(coeffs::MVector{V,Float64}, order::Int) where {V}
        #resize_coeffs1!(coeffs, order)
        return new{V}(coeffs, order)
    end
end
getindex(a::TaylorSA, n::Int64) = a.coeffs[n+1]
setindex!(a::TaylorSA, x::Float64, n::Int)  = a.coeffs[n+1] = x
function addT2(a::TaylorSA, b::TaylorSA,c::Float64,cache::TaylorSA) #where {T<:Number}
    #cache.coeffs.=a.coeffs  #not needed for in middle ops
    cache.coeffs[1]=a.coeffs[1]+b.coeffs[1]
    cache.coeffs[2]=a.coeffs[2]+b.coeffs[2]
    cache.coeffs[3]=a.coeffs[3]+b.coeffs[3]
    #cache[0]=cache[0]+c 
   # @__dot__ cache.coeffs = (+)(cache.coeffs, b.coeffs)
      return cache
end 
############################################################testing#########################################"
using BenchmarkTools
struct Data
    tSA1::TaylorSA
    tSA2::TaylorSA
    staticcache::TaylorSA
end
function createdata()
    staticCoeffs=MVector{3,Float64}(1.2,2.2,3.0)
    staticCoeffs2=MVector{3,Float64}(1.0,2.0,3.0)
    cachestaticCoeffs2=MVector{3,Float64}(0.0,0.0,0.0)
    tSA1= TaylorSA(staticCoeffs,2)
    tSA2=TaylorSA(staticCoeffs2,2)
    staticcache=TaylorSA(cachestaticCoeffs2,2)
    data=Data(tSA1,tSA2,staticcache)
end
function integrate(d::Data)
    tSA1= d.tSA1
    tSA2=d.tSA2
    staticcache=d.staticcache
    for i=1:1
    tSA1=addT2(tSA1,tSA2,1.0,staticcache)
    end
    return nothing
end
function outertest()
    d=createdata()
    integrate(d)
    return 5.0
end
#= @show testnormal()
@show teststatic() =#
#@btime testnormal()#105.558 ns (3 allocations: 240 bytes)
@btime outertest()#26.723 ns (3 allocations: 96 bytes)