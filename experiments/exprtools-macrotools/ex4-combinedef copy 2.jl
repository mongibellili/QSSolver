#module module1
using ExprTools
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
struct NLprob
    f::Expr
end
macro ode()
#=     Base.remove_linenums!(odeExpr) =#  
equs=[:(u[2]),:(u[1]+u[2]),:(-u[1]-u[3])]
@show equs
    io = IOBuffer()
    write(io, "if j==1 u[2]")
    for i=2:length(equs)-1
        write(io, "elseif j==$i $(equs[i])")
    end
    write(io, "else  $(equs[length(equs)]) end")
    s = String(take!(io))
    @show s
    close(io)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :f1
    def[:args] = [:(j::Int),:(u::Vector{Float64})]
    def[:body] = Meta.parse(s)
   
    ff=combinedef(def)
   prob=NLprob(ff)
end
ff=@ode
#@show ff.f
g=@RuntimeGeneratedFunction(ff.f)
@show g
u=[1.0,2.0,3.0]
@show g(3,u)
#= function prepareProb(eqs::Vector{Expr})
    # display(eqs[1].args[2])
    arr=Vector{Function}()
    for i=1:length(eqs)
     push!(arr,eval(eqs[i])) #f_1
    end
   arr
 end =#
#=  function integrate(u::Vector{Float64})
     
     for k=1:10
         for i=1:length(u)
             derivate(i,u,f1)
         end
     end
     #println(u)
 end
 function derivate(i::Int,u::Vector{Float64},f::Function)
     u[i]=f(u)
 end
end#end module
#-------------------------user space-----------------------------------
using BenchmarkTools
module1.@ode begin
    u[1]+2.0u[2]   
    u[2]-3
end  
u=[2.2,3.4]
#eqs=module1.prepareProb(odeprob)
@btime module1.integrate(u)
 =#
