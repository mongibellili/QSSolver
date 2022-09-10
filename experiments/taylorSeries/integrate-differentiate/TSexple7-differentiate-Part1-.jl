module modqss
    using TaylorSeries
    using StaticArrays
    using TimerOutputs
   #use_show_default(true)
   
   function qss_Integrate(initCond::Vector{Float64},f::SVector{N, Function},order::Int)where{N}
       # reset_timer!()
        states=length(initCond)
        x=Taylor1{Float64}[]
        q=Taylor1{Float64}[]
        t=Taylor1(order-1)
        for i=1:states
            push!(x,Taylor1(zeros(order+1),order)+initCond[i])
            push!(q,Taylor1(zeros(order),order-1))
        end
         for k=1:order
           # println("k= ",k)
            for i=1:states  
                q[i].coeffs[k] =x[i].coeffs[k]
            end
          # println("k $k q= ",q) 
           for i=1:states 
            x[i].coeffs[k+1] =((differentiate(f[i](q,t),k-1)).coeffs[1] )/factorial(k) 
           end
        end            
        println("initialX with derivatives= ",x)
        println("initialq with derivatives= ",q)
    end#end integrate
end#end module

#-----------user space----------------
using TaylorSeries
using StaticArrays
using BenchmarkTools
initCond=[1.0,2.4] #[x1zero,x2zero]
order=3
function f1(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    q[2]*exp(t) #cos= 1.0 - 0.5 t² + 𝒪(t³) for order 2
end
function f2(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    -q[1]-q[2]
end
jacobian=SVector{2,Function}(f1,f2) 
modqss.qss_Integrate(initCond,jacobian,order)
#@btime modqss.qss_Integrate(initCond,jacobian,order)
