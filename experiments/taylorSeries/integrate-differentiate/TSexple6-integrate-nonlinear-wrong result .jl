module modqss
    using TaylorSeries
    using StaticArrays
    using TimerOutputs
   # use_show_default(true)
   
   function qss_Integrate(initCond::Vector{Float64},f::SVector{N, Function},order::Int)where{N}
       # reset_timer!()
        states=length(initCond)
        xarr=Taylor1{Float64}[]
        t=Taylor1(order)
        for i=1:states
            push!(xarr,Taylor1(zeros(order+1),order)+initCond[i])
        end
        q=deepcopy(xarr)
         for i=1:states       
                   xarr[i] = integrate(f[i](q,t),initCond[i])        
        end       
        for k in 1:order-1  
            #@timeit "q[i].coeffs=" 
             for i=1:states
                q[i].coeffs[2:end] .= 0.0
                q[i].coeffs[1] = xarr[i].coeffs[k+1]
            end
            for i=1:states  
                xarr[i].coeffs[k+2] =integrate(f[i](q,t)).coeffs[2]
            end
        end  
        println("initialX with derivatives= ",xarr)
        #print_timer()
    end#end integrate
end#end module

#-----------user space----------------
using TaylorSeries
using StaticArrays
using BenchmarkTools

initCond=[1.0,2.4] #[x1zero,x2zero]
order=3

function f1(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    q[2]*exp(t)             #cos= 1.0 - 0.5 t² + 𝒪(t³) for order 2
end

function f2(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    -q[1]-q[2]
end

jacobian=SVector{2,Function}(f1,f2) 
modqss.qss_Integrate(initCond,jacobian,order)
#@btime modqss.qss_Integrate(initCond,jacobian,order)
