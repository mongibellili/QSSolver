module module1
using StaticArrays
#function init()

    macro ode(odeExpr)
        Base.remove_linenums!(odeExpr)   
        exs = [quote  $(Symbol(:f_, i))(u::Vector{Float64}) = $(odeExpr.args[i]) end for i = 1:2]
       # exs
       esc(Expr(:block, exs...)) 
       # (Expr(:block, exs...))
    end
    
    function prepareProb(n::Int64)
        all=true
imported=false
println(names(Main;  all,   imported))
println(Main.f_1)
       # display(eqs[1].args[2])
      # println(methods(Main.f_1))

#println(Main.module1.f_2)
     #=  arr=Vector{Function}()
       for i=1:n
        s=Symbol(:f_,i)
        #
       # const
       # f=getfield(Main,s)
        f=getfield(@__MODULE__,s)
        push!(arr,f) 
       end
       sarr=SVector{n,Function}(arr)
      sarr =#
    end
    function integrate(u::Vector{Float64},f::SVector{n,Function}) where {n}
        #display((methods(Main.f_1)))
       println((f[1]))
        for k=1:1
            for i=1:length(u)
                
                f[i](u)
            end
        end
       # println(u)
       
    end

end#end module
#-------------------------user space-----------------------------------
using BenchmarkTools
function test()
s= module1.@ode begin
    u[1]=u[1]+2.0u[2]   
    u[2]=u[2]-3
end  
u=[2.2,3.4]
#= println(methods(f_1))
println(f_2(u))  =#
all=true
imported=false
#println(names(Main;  all,   imported))
#eqs=module1.prepareProb(length(u))
#@btime 
#module1.integrate(u,eqs)
#@btime  
#println(f_1)
 #= module1.integrate(u,eqs)
display(u) =#

end
#@btime
 test()
 
