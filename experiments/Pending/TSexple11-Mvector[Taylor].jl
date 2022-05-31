module modqss
    using TaylorSeries
    using StaticArrays
    #use_show_default(true)
   
   function qss_Integrate(initCond::Vector{Float64},f::SVector{N, Function},order::Int)where{N}
        states=length(initCond)
        arr=Taylor1{Float64}[]
        for i=1:states
            xtemp=Taylor1(zeros(order+1),order)+initCond[i]
            push!(arr,xtemp)
        end
        x=MVector{states,Taylor1{Float64}}(arr) #obtain x with intiConds      
       # println("initialX without derivatives= ",x)
        q=deepcopy(x)
        #arr=Taylor1{Float64}[]        
        for i=1:states       
            xtemp = integrate(f[i](q))    
            x[i].coeffs[2] =xtemp.coeffs[2] 
            #push!(arr,xtemp)
        end       
        for k in 1:order-1  
            for i=1:states
                q[i].coeffs[2:end] .= 0.0
                q[i].coeffs[1] = x[i].coeffs[k+1]
            end
            for i=1:states  
                xtemp= integrate(f[i](q))     
                x[i].coeffs[k+2] =xtemp.coeffs[2]
            end
            #println("initialX with derivative $(k+1)= ",xarr)
        end 
       # println("initialX with derivatives= ",x) 
    end#end integrate
end#end module


#-----------user space----------------

using TaylorSeries
using StaticArrays
using BenchmarkTools
initCond=[1.0,2.4] #[x1zero,x2zero]
order=3

function f1(q::MVector{N,Taylor1{Float64}})where{N}
    q[2]+0.62
end

function f2(q::MVector{N,Taylor1{Float64}})where{N}
    -q[1]-q[2]
end
jacobian=SVector{2,Function}(f1,f2) 
#modqss.qss_Integrate(initCond,jacobian,order)
@btime modqss.qss_Integrate(initCond,jacobian,order)
