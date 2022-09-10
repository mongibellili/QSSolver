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
        return new{V}(coeffs, order)
    end
end
getindex(a::TaylorSA, n::Int64) = a.coeffs[n+1]
setindex!(a::TaylorSA, x::Float64, n::Int)  = a.coeffs[n+1] = x
function (+)(a::TaylorSA, b::Float64) 
    coeffs = copy(a.coeffs) # removing this will change also a !!!?? a.coeffs passed by reference. copy is needed, so there must be a cache.
    @inbounds coeffs[1] = (+)(a[0], b)
    return TaylorSA(coeffs, a.order)
end
############################################################testing#########################################"
using BenchmarkTools
vect=(1.2,2.2,3.0)
coeffs=[vect...]
coeffs2=[1.0,2.0,3.0]
staticCoeffs=MVector{3,Float64}(vect)
function testnormal(a::Taylor0)
    for i=1:10
         a=a+2.0
    end
    a
end
function testSA(a::TaylorSA)#where {V} 
    for i=1:10
         a=a+2.0
    end
    a
end
t1=Taylor0(coeffs,2)
tSA1= TaylorSA(staticCoeffs,2)

#= @show testnormal(t1)
@show testSA(tSA1) =#
@btime testnormal(t1) #300.211 ns (10 allocations: 800 bytes)
@btime testSA(tSA1) # 90.704 ns (11 allocations: 352 bytes)
