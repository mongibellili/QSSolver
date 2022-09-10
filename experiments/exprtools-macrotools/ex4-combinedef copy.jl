module module1
using ExprTools

macro ode(odeExpr)
    Base.remove_linenums!(odeExpr)   
    equation1=(odeExpr.args[1])
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :f1
    def[:args] = [:(u::Vector{Float64})]
    def[:body] = quote $equation1 end
    combinedef(def)
end

function f2(u::Vector{Float64})
    u[2]=u[2]-3
   end 
function integrate(u::Vector{Float64})
    display((methods(f1)));println()
   # println(f1([2.2,3.4]))
    for k=1:10
        for i=1:length(u)
            if i==1
                 f1(u)
            else
              #  f2(u)
            end
        end
    end

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
    #@btime 
     module1.integrate(u)
    display(u)
end
@btime test()

#display((methods(f1)));println()#UndefVarError: f1 not defined