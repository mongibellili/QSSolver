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
    coeffs = copy(a.coeffs) # removing this will change also a !!!?? a.coeffs passed by reference. copy is needed, so there must be a cache if removed
    @inbounds coeffs[1] = (+)(a[0], b)
    return TaylorSA(coeffs, a.order)
end
############################################################testing#########################################"
using BenchmarkTools
vect=(1.2,2.2,3.0)
coeffs=[vect...]
coeffs2=[1.0,2.0,3.0]
staticCoeffs=MVector{3,Float64}(vect)
#coeffs2=[1.2,2.2,3.0]
function testnormal(coeffs::Vector{Float64})
    t1=Taylor0(coeffs,2)
    for i=1:10
     t1=t1+2.0
    end
    t1
end
function testSA(coeffs::MVector{3,Float64})#where {V}
    tSA1= TaylorSA(coeffs,2)
    for i=1:10
     tSA1=tSA1+2.0
    end
    tSA1
end
#= @show testnormal(coeffs)
@show testSA(staticCoeffs) =#
@btime testnormal(coeffs) #303.841 ns (10 allocations: 800 bytes)
@btime testSA(staticCoeffs) # 69.030 ns (10 allocations: 320 bytes)
