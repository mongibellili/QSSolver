module module1
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
    macro ode(odeExpr)
        Base.remove_linenums!(odeExpr)   
        #exs = [quote $(Symbol(:f_, i))(u::Float64,v::Float64) = $(odeExpr.args[i]) end for i = 1:2]
        exs = [quote  $(odeExpr.args[i]) end for i = 1:2]
        exs
       # esc(Expr(:block, exs...))
    end
    
   @generated function prepareProb(eqs::Vector{Expr})
       # display(eqs[1].args[2])
      # arr=Vector{Function}()
      code = Expr(:block)
       for i=1:2
                push!(code.args, quote
                                function f(v::Val{$i})
                                eqs[i]
                                end
                                  end)
     
        end
        code
    end
#=     function integrate(x::Vector{Float64},u::Float64,v::Float64)
        #display((methods(f_1)))
        
        for k=1:1
            for i=1:2
                x[i]=f(Val(i))(u,v)
            end
        end
       # println(u)
     #  @show n
return nothing
    end =#


end#end module
#-------------------------user space-----------------------------------
using BenchmarkTools
#function test()
#= display(@macroexpand odeprob=module1.@ode begin
    u+2*v
    u-3
end )  =#

 odeprob=module1.@ode begin
    u+2*v
    u-3
end 
#x=Vector{Float64}(0.0,0.0)
x=[0.0,0.0]

u=2.2
v=3.4
module1.prepareProb(odeprob)
#@btime 
#module1.integrate(x,u,v)
#module1.integrate(u,eqs)
@show f(Val(1))
#@show x
#end
#@btime test()

