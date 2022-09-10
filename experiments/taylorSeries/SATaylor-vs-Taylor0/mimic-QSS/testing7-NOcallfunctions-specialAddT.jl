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
    #= @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs)
    cache[0]=cache[0]+c  =#
    cache.coeffs[1]=a.coeffs[1]+b.coeffs[1]+c
    cache.coeffs[2]=a.coeffs[2]+b.coeffs[2]
    cache.coeffs[3]=a.coeffs[3]+b.coeffs[3]
    return cache
end
function subT(a::Taylor0, b::Taylor0,cache::Taylor0) 
     #@__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
     cache.coeffs[1] = a.coeffs[1]-b.coeffs[1]
     cache.coeffs[2] = a.coeffs[2]-b.coeffs[2]
     cache.coeffs[3] = a.coeffs[3]-b.coeffs[3]
       return cache
 end 
################################################################SA taylor##########################################
using StaticArrays
struct TaylorSA{N} #<: AbstractSeries
    coeffs :: MVector{N,Float64}
    order :: Int
    function TaylorSA(coeffs::MVector{N,Float64}, order::Int) where {N}
        return new{N}(coeffs, order)
    end
end
getindex(a::TaylorSA, n::Int64) = a.coeffs[n+1]
setindex!(a::TaylorSA, x::Float64, n::Int)  = a.coeffs[n+1] = x
function addT2(a::TaylorSA, b::TaylorSA,c::Float64,cache::TaylorSA) #where {T<:Number}
    @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs)
    cache[0]=cache[0]+c 
    return cache
end 
function subT2(a::TaylorSA, b::TaylorSA,cache::TaylorSA) 
     @__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
       return cache
 end 
############################################################functions of TSA#########################################"
struct Data
    x::TaylorSA
    tSA::Vector{TaylorSA}
    staticcache::TaylorSA
   # farray::Vector{Function}
end
function createdata()
    staticCoeffs=MVector{3,Float64}(1.2,2.2,3.0)
    staticCoeffs2=MVector{3,Float64}(1.0,2.0,-1.35)
    cachestaticCoeffs2=MVector{3,Float64}(0.0,0.0,0.0)
    xCoeffs2=MVector{3,Float64}(0.0,0.0,0.0)
    x=TaylorSA(xCoeffs2,2)
    tSA= [TaylorSA(staticCoeffs,2),TaylorSA(staticCoeffs2,2)]
    #tSA2=TaylorSA(staticCoeffs2,2)
    staticcache=TaylorSA(cachestaticCoeffs2,2)
    data=Data(x,tSA,staticcache)
end
 function f1(x::TaylorSA,j::Int,tSA::Vector{TaylorSA},staticcache::TaylorSA)#@inline hurts performance
    x
    if j==1
    x=addT2(tSA[1],tSA[2],1.0,staticcache)
   # @show x
    else
    x=subT2(tSA[2],tSA[1],staticcache)
    end
    return nothing
end
function integrate(d::Data)
    x=d.x
    #f=d.farray
    tSA= d.tSA
    staticcache=d.staticcache
    index=1
    for i=1:1
        f1(x,index,tSA,staticcache)
        index=3-index
    end
    #return nothing
    return x
end
function outertest()
    d=createdata()
    integrate(d)
    #display(@code_typed integrate(d))
end
############################################################functions of Taylor0#########################################"
using BenchmarkTools
struct Data0
    x::Vector{Taylor0}
    q::Vector{Taylor0}
    cache::Taylor0
end
function createdata0()
    coeffs=[1.2,2.2,3.0]
    coeffs2=[1.0,2.0,-1.35]
    cacheCoeffs=[0.0,0.0,0.0]
    x1coeffs=[0.0,0.0,0.0]
    x2coeffs=[0.0,0.0,0.0]
    x=[Taylor0(x1coeffs,2),Taylor0(x2coeffs,2)]
    q=[Taylor0(coeffs,2),Taylor0(coeffs2,2)]
    cache=Taylor0(cacheCoeffs,2)
    data=Data0(x,q,cache)
end
#= function f2(x::Vector{Taylor0},j::Int,q::Vector{Taylor0},cache::Taylor0)
    if j==1
        temp=addT(q[1],q[2],1.0,cache)
        x[1][0]=temp[0]
        x[1][1]=temp[1]
        x[1][2]=temp[2]

    else
       temp=subT(q[2],q[1],cache)
       x[2][0]=temp[0]
       x[2][1]=temp[1]
       x[2][2]=temp[2]


    end
    return nothing
end =#
function integrate0(d::Data0)
    x=d.x
    q=d.q
    cache=d.cache
    index=1
    for i=1:100
        if index==1
            temp=addT(q[1],q[2],1.0,cache)
            x[1][0]=temp[0]
            x[1][1]=temp[1]
            x[1][2]=temp[2]
    
        else
           temp=subT(q[2],q[1],cache)
           x[2][0]=temp[0]
           x[2][1]=temp[1]
           x[2][2]=temp[2]
    
    
        end
        index=3-index
    end
    #return nothing
    return x
end
function outertest0()
    d=createdata0()
    integrate0(d)
end

###############################################testing##############################""
using BenchmarkTools
using InteractiveUtils

#@show outertest()
#@show outertest0()

#= outertest() = TaylorSA([-0.19999999999999996, -0.20000000000000018, -4.35], 2)
outertest0() = Taylor0([-0.19999999999999996, -0.20000000000000018, -4.35], 2) =#


#@btime outertest()# 1.843 μs (30 allocations: 1.06 KiB)
@btime outertest0()#1.738 μs (26 allocations: 1.09 KiB)
