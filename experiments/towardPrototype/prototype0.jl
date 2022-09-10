

 module module1
 using StaticArrays
 using TimerOutputs
 using InteractiveUtils
 using RuntimeGeneratedFunctions
 RuntimeGeneratedFunctions.init(@__MODULE__)

  macro ode(odeExpr)
         Base.remove_linenums!(odeExpr)   
         exs = [quote $(Symbol(:f_, i))(u::Vector{Float64}) =$(odeExpr.args[i]) end for i = 1:2]
         exs
  end
     

   function integrate(u::Vector{Float64},eqs::Vector{Expr})
   # reset_timer!()
    arr=Vector{Function}()
    for i=1:length(eqs)
    # push!(arr,eval(eqs[i])) #f_1
     push!(arr, @RuntimeGeneratedFunction(eqs[i].args[2]))
    end
         for k=1:1
             for i=1:length(u)
              #  @timeit "derivate"  
              arr[i](u)
             # derivate(i,u,arr[i])
              #  display(@code_typed derivate(i,u,arr[i]))
             end
         end
        # println(u)
       #print_timer() 
   end
   function derivate(i::Int,u::Vector{Float64},f::Function)
         # u[i]=
          f(u)
         return nothing
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
 
 #@btime module1.integrate(u,eqs)
 module1.integrate(u,odeprob)
 #display(u)
 end
# @btime 
 #test()
 
 #= display(@macroexpand module1.@ode begin
    u[1]=u[1]+2.0u[2]   
    u[2]=u[2]-3
end  ) =#
 


