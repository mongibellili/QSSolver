
module smallqss
using StaticArrays

#= function equations(du::MVector{R,Float64},u::MVector{R,Float64},p::Float64,t::Float64)
    Base.remove_linenums!(schema)
    #rhs=:($schema).args[1].args[2]
    v=[]
    for i=1:length(schema.args) 
        code2=(quote  
                        j=$i
                        #@inline 
                        function f(::Val{j},u::MVector{R,Float64}) where {R}
                         $(:($schema).args[:($i)].args[2])# esc(:($schema).args[1].args[2])
                       end
               end)
        push!(v,code2)      
    end
       esc(Expr(:block,v...))   
end =#

function solve(fun,u::MVector{R,Float64}) where {R} # test the created function here
    for i=1:50
        for index=1:length(u)
            computeDerivative(fun,index,u)
        end
    end
end 

function computeDerivative(fun,index::Int,u::MVector{R,Float64}) where {R} # test the created function here   
       u[index]= Main.f(Val(index),u)
      # u[index]=u[2] 
      # u[index]= ff(Val(index),u)
end 
function ODEProblem(func)
    dump(func)
end
end #end module
#--------------------------user space--------------------
using StaticArrays
using BenchmarkTools
function test()
    function func(du,u,p,t)
        du[1] = u[2] 
        du[2] = -u[1] - u[2]  
    end
    prob = smallqss.ODEProblem(func)
    display(prob);println() 

    u=@MVector[1.1,2.6]
    #= display(f(Val(1),u));println() #2.6
    display(f(Val(2),u));println() #-3.7 =#
    #prob=problem(odefun)
    #display(prob)
    #smallqss.solve() #UndefVarError: f not defined....f lives in the user world....
   # smallqss.solve(prob,u)
   # display(u)
  # smallqss.solve(u)
    #display(methods(f))
end
#@btime test()
test()
