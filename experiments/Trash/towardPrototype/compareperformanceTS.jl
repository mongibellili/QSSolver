using TaylorSeries
using BenchmarkTools
using InteractiveUtils

#############################################################################################
#                                                                                           #
#       Warning: differentiate! and evaluateX  are not exported in the original TS          #
#                                                                                           #
#############################################################################################



function differentiate!(res::Taylor1{Float64}, a::Taylor1{Float64})
    for k in eachindex(res)
       # differentiate!(res, a, ord)
        if k < a.order
          @inbounds res[k] = (k+1)*a[k+1]
      end
    end
    nothing
  end
function code1(x::Vector{Taylor1{Float64}},t::Taylor1{Float64})
#x[1]=x[1](t+0.5)
#evaluateX(x[1],t+0.5)# did not remove allocations
return nothing
end

function code2(x::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    order=3
    for k = 1:order #evaluate x at new time, and derivatives...remember coeffs store the over facto also              
      x[1].coeffs[k] = (((differentiate(x[1], k - 1)(t + 0.5))).coeffs[1]) / factorial(k - 1)
    end
    return nothing
end
function code3(x::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    order=3
    for k = 1:order            
      x[1].coeffs[k] = (((differentiate(x[1], k - 1)(0.5)))) / factorial(k - 1)
    end
   return nothing
end
function code33(x::Vector{Taylor1{Float64}},cacheT::Taylor1{Float64})
    x[1].coeffs[1] = x[1](0.5)
    differentiate!(cacheT,x[1])
    x[1].coeffs[2] = cacheT(0.5)
    #order=3
    for k = 3:3            
      differentiate!(cacheT,cacheT)
      x[1].coeffs[k] = cacheT(0.5)/2
    end
  # return nothing
end
function code4(x::Vector{Taylor1{Float64}})               
      x[1].coeffs[1] = x[1](0.5)
      #c=differentiate(x[1])
     # x[1].coeffs[2] = c(0.5)
     x[1].coeffs[2] = differentiate(x[1])(0.5)
    #  d=differentiate(c)
    #  x[1].coeffs[3] = d(0.5)/2
     # return nothing
end
function code44(x::Vector{Taylor1{Float64}},cacheT::Taylor1{Float64})     
    x[1].coeffs[1] = x[1](0.5)
    differentiate!(cacheT,x[1])
    x[1].coeffs[2] = cacheT(0.5)
   # differentiate!(cacheT,cacheT)
   # x[1].coeffs[3] = cacheT(0.5)/2
    #return nothing
end
t=Taylor1(3)
a=Taylor1([1.0,2.0,5.0,-4.4])
b= Taylor1([2.2,2.8,3.0,-2.5])
cacheT=Taylor1([0.0,0.0,0.0,0.0])
x=[a,b]
#code1(x,t)
#code2(x,t)
#code3(x,t)
#code4(x)
#code33(x,cacheT)
#code44(x,cacheT)
 # println(" :x1 after evaluate = ", x[1])
   #= @btime code1(x,t)
@btime code2(x,t) 
@btime code3(x,t) 
@btime code4(x) =#

#display(@code_lowered code4(x))
#display(@code_llvm code44(x,cacheT))
#display(@code_native code4(x))
display(@code_typed code4(x))


#= 
@btime code33(x,cacheT) 
@btime code44(x,cacheT)  =#
#display(@code_typed code1(x,t))
#display(@code_warntype code1(x,t))
#display(@code_typed code2(x,t))
#display(@code_warntype code2(x,t))
#display(@code_typed code3(x,t))
#display(@code_warntype code3(x,t))
#display(@code_typed code4(x))
#display(@code_warntype code4(x))

#


#display(@code_warntype code33(x,cacheT))

