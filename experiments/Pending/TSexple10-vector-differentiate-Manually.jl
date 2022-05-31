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
            for i=1:states  
                q[i].coeffs[1] =x[i].coeffs[1]
            end
            x[1].coeffs[2] =(f[1](q,t)).coeffs[1] 
            x[2].coeffs[2] =(f[2](q,t)).coeffs[1] 
            for i=1:states  
                q[i].coeffs[2] =x[i].coeffs[2]
            end
            x[1].coeffs[3] =(Main.derf1(q,t).coeffs[1] )/factorial(2) 
            x[2].coeffs[3] =(differentiate(f[2](q,t),1).coeffs[1])/factorial(2)           
            for i=1:states  
                q[i].coeffs[3] =x[i].coeffs[3]
            end
            println("q= ",q) 
            x[1].coeffs[4] =(Main.derderf1(q,t).coeffs[1] )/factorial(3) 
            x[2].coeffs[4] =(differentiate(f[2](q,t),2).coeffs[1]) /factorial(3)
            #println("initialX with derivative $(k+1)= ",xarr)      
        println("initialX with derivatives= ",x)
       # println("initialq with derivatives= ",q)
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
    q[2]*exp(t) #cos= 1.0 - 0.5 t² + 𝒪(t³) for order 2
end
function derf1(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    differentiate(q[2])*exp(t)+ q[2]*differentiate(exp(t))
end
function derderf1(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    a1=(differentiate(q[2],2))#/factorial(3)
    println("a1= ",a1)
    a2=exp(t)
    println("a2= ",a2)
    a=a1*a2
    println("a= ",a)
    b1= differentiate(q[2])
    println("b1= ",b1)
    b2=differentiate(exp(t))
    println("b2= ",b2)
    b=b1*b2
    println("b= ",b)
    c1=q[2]
    println("c1= ",c1)
    c2=differentiate(exp(t),2)#/factorial(3)
    println("c2= ",c2)
    c=c1*c2
    println(c)
    a+2*b+c

end
function f2(q::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    -q[1]-q[2]
end
jacobian=SVector{2,Function}(f1,f2) 
modqss.qss_Integrate(initCond,jacobian,order)
#@btime modqss.qss_Integrate(initCond,jacobian,order)
