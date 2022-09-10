using TaylorSeries
using BenchmarkTools
a =   [ Taylor1([1.0,1.1,1.2],2),Taylor1([1.3,1.4,1.5],2)]
b =   [ Taylor1([2.0,2.1,2.2],2),Taylor1([2.3,2.4,2.5],2)]
c =   [ Taylor1([3.0,3.1,3.2],2),Taylor1([3.3,3.4,3.5],2)]
d =    [Taylor1([4.0,4.1,4.2],2),Taylor1([4.3,4.4,4.5],2)]
cache=Taylor1([0.0,0.0,0.0],2)

function normalArith(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cacheT::Taylor1{Float64})
    #sum=a[1]-b[1]-c[1]-d[1]+a[2]-b[2]
    sum= b[1]- a[1]
    return sum
  end
   @btime normalArith(a,b,c,d,cache)

 #=   function stressTestMul1(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cache::Vector{Taylor1{Float64}})
    #test 2 vars
    @assert a[1]*1.3==@changeAST a[1]*1.3,12
    @show cache[1]
    @show cache[2]
    @assert 1.3*a[1]==@changeAST 1.3*a[1],12
    @show cache[1]
    @show cache[2]
    @assert a[1]*b[1]==@changeAST a[1]*b[1],12
    @show cache[1]
    @show cache[2]
    @assert 3.2*1.3==(@changeAST 3.2*1.3,12)[0]
    @show cache[1]
    @show cache[2]
    @assert a[1]*1.3*2.0==@changeAST a[1]*1.3*2.0,12
    @show cache[1]
    @show cache[2]
    @assert 1.3*a[1]*2.0==@changeAST 1.3*a[1]*2.0,12
    @show cache[1]
    @show cache[2]
    @assert a[1]*b[1]*2.0==@changeAST a[1]*b[1]*2.0,12
    @show cache[1]
    @show cache[2]
     #test 3 vars
     @assert a[1]*1.3*1.2==@changeAST a[1]*1.3*1.2,12
     @show cache[1]
    @show cache[2]
     @assert 1.3*a[1]*1.2==@changeAST 1.3*a[1]*1.2,12
     @show cache[1]
    @show cache[2]
     clearCache(cache)
     @assert 1.3*1.2*a[1]==@changeAST 1.3*1.2*a[1],12
     @show cache[1]
    @show cache[2]
     @assert a[1]*b[1]*1.2==@changeAST a[1]*b[1]*1.2,12
     @show cache[1]
    @show cache[2]
     @assert a[1]*1.2*b[1]==@changeAST a[1]*1.2*b[1],12
     @show cache[1]
    @show cache[2]
     @assert a[1]*b[1]*1.2==@changeAST a[1]*b[1]*1.2,12
     @show cache[1]
     @show cache[2]
     @assert a[1]*b[1]*c[1]==@changeAST a[1]*b[1]*c[1],12
     @show cache[1]
    @show cache[2]
     @assert 1.3*1.2*2.3==(@changeAST 1.3*1.2*2.3,12)[0]
    
 end =#
 