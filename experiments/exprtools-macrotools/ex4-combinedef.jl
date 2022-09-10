module module1
#using ExprTools
#= macro ode(odeExpr)
    Base.remove_linenums!(odeExpr)   
    equation1=(odeExpr.args[1])
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :f1
    def[:args] = [:(u::Vector{Float64})]
    def[:body] = quote $equation1 end
    combinedef(def)
end =#

#= function prepareProb(eqs::Vector{Expr})
    # display(eqs[1].args[2])
    arr=Vector{Function}()
    for i=1:length(eqs)
     push!(arr,eval(eqs[i])) #f_1
    end
   arr
 end =#
 function integrate(u::Vector{Float64},f::Function)
     du=[0.0,0.0]
     for k=1:1000
         for i=1:length(u)
             #derivate(i,du,u,f)
            f(du,u)
         end
     end
     #println(du)
 end
 @inline  function derivate(i,du,u,f)
     f(du,u)
    
 end
end#end module
#-------------------------user space-----------------------------------
using BenchmarkTools
function ode(du,u)
    du[1]=u[1]+2.0u[2]   
    du[2]=u[2]-3
    return nothing
end  
u=[2.2,3.4]
#eqs=module1.prepareProb(odeprob)
@btime module1.integrate(u,ode)
#module1.integrate(u,ode)

