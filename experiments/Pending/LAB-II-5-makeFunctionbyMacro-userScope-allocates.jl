
module smallqss
using StaticArrays

macro equations(schema)
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
end

#= @generated function Odefun(schema...)
  println(schema)
end

function problem(schema::Expr) #this calls the macro that will create the new function
end 

function solve() # test the created function here
    u=@MVector[1.1,2.6,3.1]
    #= v=f(Val(1),u)
    display(v) =#
    display(f(Val(1),u));println()
end  =#
function ff(::Val{1},u::MVector{R,Float64}) where {R}
    u[2] 
end
function ff(::Val{2},u::MVector{R,Float64}) where {R}
    -u[1]-u[2] 
end
function ff(::Val{3},u::MVector{R,Float64}) where {R}
    u[3]+u[2]
end
function solve(f::Function,u::MVector{R,Float64}) where {R} # test the created function here
    for i=1:50
        for index=1:length(u)
            computeDerivative(index,f,u)
        end
    end
end 


function computeDerivative(index::Int,f::Function,u::MVector{R,Float64}) where {R} # test the created function here   
       u[index]= f(Val(index),u)
      # u[index]=u[2] 
       #u[index]= ff(Val(index),u)
end 


end #end module
#--------------------------user space--------------------
using StaticArrays
using BenchmarkTools
myExp=smallqss.@equations begin 
    du[1]=u[2] 
    du[2]=-u[1]-u[2] 
    du[3]=u[3]+u[2]
end 
display(myExp);println() #f (generic function with 1 method)

u=@MVector[1.1,2.6,3.1]
#= display(f(Val(1),u));println() #2.6
display(f(Val(2),u));println() #-3.7 =#
#prob=problem(odefun)
#display(prob)
#smallqss.solve() #UndefVarError: f not defined....f lives in the user world....
#smallqss.solve(myExp,u)
display(u)
@btime smallqss.solve(myExp,u)
#display(methods(f))