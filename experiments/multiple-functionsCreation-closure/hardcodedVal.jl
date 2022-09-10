module module1
    function integrate(u::Vector{Float64})      
        for k=1:10
            for i=1:length(u)             
                f(Val(i),u)
            end
        end
       # println(u)
       
    end
function f(::Val{1},u::Vector{Float64})
    u[1]=u[1]+2.0u[2] 
end
function f(::Val{2},u::Vector{Float64})
    u[2]=u[2]-3 
end
end#end module
#-------------------------user space-----------------------------------
using BenchmarkTools
#function test()
#= module1.@ode begin
    u[1]=u[1]+2.0u[2]   
    u[2]=u[2]-3
end  =# 
u=[2.2,3.4]

@btime module1.integrate(u)
#display(u)

#end
#@btime
 #test()

