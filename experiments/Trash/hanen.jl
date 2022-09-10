#= 
 test(a, b, c, xs...)=test2(a,b,xs...)

function test2(a,c,d...)
  @show length(d)
 @show d
 @show a
 
 @show c
end
test(1,2,3)
 =#

#z=:(sqrt(5))#z=:(exp(5))#z=:(log(5))#z=:(abs(5))
#= head: Symbol call
  args: Array{Any}((2,))
    1: Symbol sqrt
    2: Int64 5 =#
#z=:(5^2)
#= head: Symbol call
  args: Array{Any}((3,))
    1: Symbol ^
    2: Int64 5
    3: Int64 2 =#
#= z=:(abs(5))
dump(z) =#
#= z=:(p=5.0)
dump(z) =#
#= y=:(z+0.0)
dump(y) =#
#= cachexpr1 = Expr(:call, :+,z,0.0)
dump(cachexpr1) =#

#= 
using InteractiveUtils
using BenchmarkTools

function addsub(a::String,c::Int,b::String)
  # a*sting(c)*b
    string(a,c,b)
end
function addsub(a::Int,c::String,b::String)
      string(a,c,b)
  end
  function addsub(a::String,c::String,b::Int)
    # a*sting(c)*b
      string(a,c,b)
  end
  function addsub(a::Int,c::Int,b::String)
    # a*sting(c)*b
      string(a,c,b)
  end
  function addsub(a::String,c::Int,b::Int)
    # a*sting(c)*b
      string(a,c,b)
  end
  function addsub(a::Int,c::String,b::Int)
    # a*sting(c)*b
      string(a,c,b)
  end

  function addsub(a::Int,c::Int,b::Int)
    # a*sting(c)*b
      string(a,c,b)
  end
  function addsub(a::String,c::String,b::String)
    # a*sting(c)*b
      string(a,c,b)
  end

function subadd(a::P,b::Q,c::R) where {P,Q,R <:Union{String,Int}}
    addsub(a::P,c::R,b::Q)  
end



#println(@code_warntype subadd("1",2,3))
#= @show subadd("1","2",3)
@show subadd(1,"2",3)
@show subadd("1",2,3) =#
@btime subadd("1",2,3) =#

#= using TaylorSeries
#use_show_default(true)
t0=Taylor1([4.0,2.0,0.0],2) =#
#t0=Taylor1([0.0im,1.0,0.0],2)
#= t1=Taylor1([0.0,1.0,0.0],2)
t2=Taylor1([1.0,1.0,1.0],2) =#
#@show t0
#@show log(t1)
#@show sincos(3.14)
#@show abs(t0)
#@show sqrt(abs2(t0))
#m=abs2(t0)
#n=Taylor1{Float64}([0.0, 1.0, 1.0], 2)
#@show abs2(t0)
#@show sqrt(n)
#@show sqrt(t0)
#@show sqrt(t1)
#p=1
#= for p=1:5
 #println("for p= ",p) 
t = trailing_zeros(p) +1

#@show t
        @show p >> t
#@show p
end =#
#@show t0^0.49999999999
#@show sqrt(t0)
#= c1 = Taylor1( 0.0, 2 )
c1[0]=sqrt(4.0)
@show c1
c2 = Taylor1( sqrt(4.0), 2 )
@show c2
@show sqrt(t0) =#
#= 
a=1.1213
b=1.1213
r=2*(2*a+b)^3-9*((2*a+b)*(2*a*b+a*a))+27*a*a*b
@show r =#
#= using BenchmarkTools
function f(v::SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true})
  #display(v) 
  res=v[1]
end
function f(v::Vector{Float64})
  #display(v) 
  res=v[1]
end
function f(v::Base.ReshapedArray{Float64, 1, Vector{Float64}, Tuple{}})
  #display(v) 
  res=v[1]
end
N=5
a2 = rand(N)
@show a2
#= @show f(@view a2[1:3])
@btime f(@view a2[1:3]) =#
@show f(  reshape(a2[2:3],2)  )
@show f(  Base.ReshapedArray((a2),(1,),())  )
@btime f(  Base.ReshapedArray((a2),(2,),())  ) =#
#= 
function f(x)
  cos(x)^3-0.75cos(x)-0.25cos(3x)
end
@show f(0.225)

function f2(x)
  6cos(x)+2cos(3x)
end
@show f2(0.225) =#

#= result = Ref{Int64}(0)
t1=ccall((:time, "libc.so.6"), Int64, (Ptr{Int64},), result)
#@show t1 #1662311444
#@show result#Base.RefValue{Int64}(1662311425) #result[]=1662311444  
mutable struct Ctm
  sec::Cint
  min::Cint
  hour::Cint
  mday::Cint
  mon::Cint
  year::Cint
  wday::Cint
  yday::Cint
  isdst::Cint
end
function Base.show(io::IO, t::Ctm)
  print(io, t.year + 1900, "-", lpad(t.mon, 2, "0"), "-",
        lpad(t.mday, 2, "0"), "T", lpad(t.hour, 2, "0"), ":",
        lpad(t.min, 2, "0"), ":", lpad(t.sec, 2, "0"))
end
localtime = ccall((:localtime, "libc.so.6"), Ptr{Ctm}, (Ptr{Int64},), result)
#@show localtime  #Ptr{Ctm} @0x00007f7adc66a780
t2=unsafe_load(localtime)
#@show t2  #2022-08-04T19:16:25 =#

mutable struct Bar
  x::Int
end
b = Bar(2)
#@show b
myptr=pointer_from_objref(b)
#@show myptr #Ptr{Nothing} @0x00007f9234db4b10
myptr1 = Ptr{Bar}()#Ptr{Bar} @0x0000000000000000
myptr2 = Ptr{Bar}(myptr)#@0x00007fec2fdf0b10
@show myptr2
t2=unsafe_load(myptr2) #Bar(2)
@show t2
