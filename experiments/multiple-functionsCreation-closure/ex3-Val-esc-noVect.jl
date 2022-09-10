module module1
#function init()

    macro ode(odeExpr)
        Base.remove_linenums!(odeExpr)   
        exs = [quote $(Symbol(:f_, 1))(::Val{$i},u::Vector{Float64}) = $(odeExpr.args[i]) end for i = 1:2]
       # exs
        esc(Expr(:block, exs...))
    end
    
   #=  function prepareProb(eqs::typeof(f_2))
       # display(eqs[1].args[2])
       arr=Vector{Function}()
       for i=1:length(eqs)
        push!(arr,eval(eqs[i])) #f_1 # this does not work when user puts code inside another function
       # push!(arr, @RuntimeGeneratedFunction(eqs[i].args[2]))
       end
      arr
    end =#
    function integrate(u::Vector{Float64})#,f::Vector{Function})
        #display((methods(Main.f_1)))
       
        for k=1:1
            for i=1:length(u)
              
                Main.f_1(Val(i),u)
            end
        end
       # println(u)
       
    end

end#end module
#-------------------------user space-----------------------------------
using BenchmarkTools
function test()
module1.@ode begin
    u[1]=u[1]+2.0u[2]   
    u[2]=u[2]-3
end  
u=[2.2,3.4]

 module1.integrate(u)
display(u)

end
#@btime
 test()

