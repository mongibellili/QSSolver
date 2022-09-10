#= 
syms = [gensym() for _ in 1:5]
syms[1] = Symbol(:1.0)
@show syms =#
#= f(x)=x*x-2x-1
@show f(2.414213562506035) =#
#= using Plots
#gr();
f(x) = x*x*x-1.5x*x-0.25x-2.5
derf(x) = 3*x*x-3x-0.25


display(plot!(f,label="pol",color=:orange))
display(plot!(derf,label = "pol'",color = :purple, xlims = (0.0,4.0),ylims = (-4.0,10.0)))
println("done")
        readline() =#

        @show typemax(Float64)
        @show  typemin(Float64)