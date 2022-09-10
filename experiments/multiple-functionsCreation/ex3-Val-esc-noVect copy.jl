module module1
#function init()

    macro ode(odeExpr)
        Base.remove_linenums!(odeExpr)  
        T=2 
        code = Expr(:block)
        for i in 1:T
            if i==1
                codepart=quote   $(Symbol(:f_, 1))(::Val{1},u::Vector{Float64}) = $(odeExpr.args[i]) end
            elseif i==2
                 codepart=quote $(Symbol(:f_, 1))(::Val{2},u::Vector{Float64}) = $(odeExpr.args[i]) end
            else
                codepart=quote     $(Symbol(:f_, 1))(::Val{3},u::Vector{Float64}) = $(odeExpr.args[i]) end
               
            end
          push!(code.args, codepart)
        end
        println(code)
        esc(code)
      #=   exs = [quote $(Symbol(:f_, 1))(::Val{$i},u::Vector{Float64}) = $(odeExpr.args[i]) end for i = 1:T]
       # exs
        esc(Expr(:block, exs...)) =#
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
       
        for k=1:10
            for i=1:length(u)
              
                Main.f_1(Val(i),u)
            end
        end
       # println(u)
       
    end

end#end module
#-------------------------user space-----------------------------------
using BenchmarkTools
#function test()
module1.@ode begin
    u[1]=u[1]+2.0u[2]   
    u[2]=u[2]-3
end  
u=[2.2,3.4]

@btime  module1.integrate(u)
#display(u)

#end
#@btime test()

