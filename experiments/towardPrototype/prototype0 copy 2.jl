

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
    reset_timer!()
    arr=Vector{Function}()
    #= for i=1:length(eqs)
    # push!(arr,eval(eqs[i])) #f_1
     push!(arr, @RuntimeGeneratedFunction(eqs[i].args[2]))
    end =#
    #= for i=1:length(eqs)
        println("(eqs[i].args[2])= ",(eqs[i].args[2]))
         push!(arr,eval(eqs[i].args[2])) #f_1:running in world age problem
       #  push!(arr, @RuntimeGeneratedFunction(eqs[i].args[2]))
        end =#
    #= push!(arr,f1)
    push!(arr,f2) =#

         for k=1:10
             for i=1:length(u)
              #  @timeit "derivate"  
             @timeit  "arr[i]" arr[i](u)  #3.81MiB     616KiB if push!(arr,f1)
          #conclusion1: @RuntimeGeneratedFunction is responsible for 3.81MiB -616KiB
            #conclusion2: arr[i] is responsible for 616KiB
             #= @timeit "derivate"   if i==1
                #f1(u) #0.00B
                arr[1](u)#3.81MiB    616KiB if push!(arr,f1)
             else
               # f2(u)
               arr[2](u)
             end  =#
             # derivate(i,u,arr[i])
              #  display(@code_typed derivate(i,u,arr[i]))
             end
         end
        # println(u)
       print_timer() 
   end
   function f1(u::Vector{Float64})
    u[1]=u[1]+2.0u[2] 
   end
   function f2(u::Vector{Float64})
    u[2]=u[2]-3
   end

  #=  function derivate(i::Int,u::Vector{Float64},f::Function)
         # u[i]=
          f(u)
         return nothing
     end =#
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
 display(u)
 end
# @btime 
 test()
 
 #= display(@macroexpand module1.@ode begin
    u[1]=u[1]+2.0u[2]   
    u[2]=u[2]-3
end  ) =#
 


