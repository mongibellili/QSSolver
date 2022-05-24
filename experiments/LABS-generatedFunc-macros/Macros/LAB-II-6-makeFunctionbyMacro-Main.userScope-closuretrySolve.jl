
module smallqss
using StaticArrays

macro equations(schema)
    Base.remove_linenums!(schema)
    #rhs=:($schema).args[1].args[2]
    v=[]
    j=0
    for i=1:length(schema.args) 
        f= Symbol(:h_,i)
        code2=(quote  
                        #j=$i
                        #@inline 
                        #scope=@__MODULE__ #__module__
                       #=  sy=Symbol("f",$i)
                        println(sy) =#
                        #f= Symbol(:h_,$i)
                        #display(@eval $f);println()
                        function $(esc(f))(u::MVector{R,Float64}) where {R}
                          #$(:($schema).args[:($i)].args[2])
                          $(:($(schema)).args[:($i)].args[2])         
                        end           
               end)
        push!(v,code2)      
    end
       (Expr(:block,v...))   
end

function solve(u::MVector{R,Float64}) where {R} # test the created function here
   # for i=1:50
        for index=1:length(u)
            computeDerivative(index,u)
        end
    #end
end 

function computeDerivative(index::Int,u::MVector{R,Float64}) where {R} # test the created function here   
       s=Symbol(:h_,index)
       #@__MODULE__
       f=getfield(Main,s)
       u[index]= f(u)
      # u[index]= f(Val(index),u)
      # u[index]=u[2] 
      # u[index]= ff(Val(index),u)
end 


end #end module
#--------------------------user space--------------------
using StaticArrays
using BenchmarkTools
function test()  # error local j cannot be used in closure ...val{j}
#display(@macroexpand
 smallqss.@equations begin 
        du[1]=u[2] 
        du[2]=-u[1]-u[2] 
        du[3]=u[3]+u[2]
    end 
   # display(myExp);println() #f (generic function with 1 method)

    u=@MVector[1.1,2.6,3.1]
    #= display(f(Val(1),u));println() #2.6
    display(f(Val(2),u));println() #-3.7 =#
    #prob=problem(odefun)
    #display(prob)
    #smallqss.solve() #UndefVarError: f not defined....f lives in the user world....
    smallqss.solve(u)
    display(u)
  # @btime smallqss.solve(u)
    display(methods(h_1))
    display(methods(h_2))
    display(methods(h_3))
end
#@btime 
test()

