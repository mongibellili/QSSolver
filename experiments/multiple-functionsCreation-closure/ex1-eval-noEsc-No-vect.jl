module module1

    macro ode(odeExpr)
        Base.remove_linenums!(odeExpr)   
        exs = [quote $(Symbol(:f_, i))(u::Vector{Float64}) = $(odeExpr.args[i]) end for i = 1:2]
        exs
       # esc(Expr(:block, exs...))
    end
    
    
    function integrate(u::Vector{Float64},eqs::Vector{Expr})
        #display((methods(f_1)))
        for k=1:1
            for i=1:length(u)
                eval(eqs[i])(u)
            end
        end
       # println(u)
       
    end

end#end module
#-------------------------user space-----------------------------------
using BenchmarkTools
function test()
odeprob=module1.@ode begin
    u[1]=u[1]+2.0u[2]   
    u[2]=u[2]-3
end  
u=[2.2,3.4]

#@btime module1.integrate(u,odeprob)

module1.integrate(u,odeprob)
display(u)

end
#@btime 
test()

