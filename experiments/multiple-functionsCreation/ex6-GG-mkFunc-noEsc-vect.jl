module module1
using GeneralizedGenerated
    macro ode(odeExpr)
        Base.remove_linenums!(odeExpr)   
        exs = [quote (u::Vector{Float64}) -> $(odeExpr.args[i]) end for i = 1:2]
        
        exs
       # esc(Expr(:block, exs...))
    end
    
    function prepareProb(eqs::Vector{Expr})
       # display(eqs[1].args[2])
       arr=Vector{Function}()
   
       println(eqs)
       println(eqs[1].args[2])
       for i=1:length(eqs)
        push!(arr,  mk_function(eqs[i].args[2]))#cannot  convert gg to function
       end
       println(arr)
      arr
    end
    function integrate(u::Vector{Float64},f::Vector{Function})
        #display((methods(f_1)))
        for k=1:1
            for i=1:length(u)
                f[i](u)
            end
        end
       # println(u)
       
    end

end#end module
#-------------------------user space-----------------------------------
using BenchmarkTools
#function test()
odeprob=module1.@ode begin
    u[1]=u[1]+2.0u[2]   
    u[2]=u[2]-3
end  
u=[2.2,3.4]
eqs=module1.prepareProb(odeprob)
#@btime module1.integrate(u,eqs)
#module1.integrate(u,eqs)
#display(u)

#end
#@btime test()

