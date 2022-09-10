module module1

    macro ode(odeExpr)
        Base.remove_linenums!(odeExpr)   
        exs = [quote function f(v::Val{$i},u::Vector{Float64}) ; $(odeExpr.args[i]) end end for i = 1:2]
       exs
        #esc(Expr(:block, exs...))
    end
    
    function prepareProb(eqs::Vector{Expr})
      
       for i=1:length(eqs)
        eval(eqs[i])
    
      # @RuntimeGeneratedFunction(eqs[i].args[2]) #yields function with not escaped unique name
     
       end
      
    end
    function integrate(u::Vector{Float64})
        #display((methods(f_1)))
        for k=1:10
            for i=1:length(u)
                f(Val(i),u)
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
module1.prepareProb(odeprob)
module1.integrate(u)
#module1.integrate(u,eqs)
#display(u)

end
#@btime 
test()

