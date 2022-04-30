using BenchmarkTools
using StaticArrays
using qss
function qssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0]
    jacobian=@SMatrix[0.0 1.0;-1.0 -1.0]
    settings = ModelSettings(initConditions,jacobian,5.5,saveat(0.12),qss1())#do not call saveat to not save, i should fix when called with zero; it does not save at all
    sol=QSS_simGenerate(settings)
end
@btime qssApproachInitialInside()
#qssApproachInitialInside()
#display(@benchmark qssApproachInitialInside())
#test plot
#= using Plots;gr()
using BenchmarkTools
using StaticArrays
using qss
function qssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0]
    jacobian=@SMatrix[0.0 1.0;-1.0 -1.0]
    settings = ModelSettings(initConditions,jacobian,10.5,saveat(0.1),qss1())#do not call saveat to not save, i should fix when called with zero; it does not save at all
    sol=QSS_simGenerate(settings)
    display(plot!(sol[1],sol[2][1]))
    display(plot!(sol[1],sol[2][2]))
end
qssApproachInitialInside()
println("done") 
readline() =#





#= #test plot with other package
using Plots;gr()
using BenchmarkTools
using StaticArrays
using qss
using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2]
        du[2] = -u[1] - u[2]   
    end
    tspan = (0.0,10)
    u0 = [1.0,2.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
    display(plot!(sol,line=(:dot, 4)))
end

function qssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0]
    jacobian=@SMatrix[0.0 1.0;-1.0 -1.0]
    #states=2
    settings = ModelSettings(initConditions,jacobian,10.5,saveat(0.1),qss1())#do not call saveat to not save, i should fix when called with zero; it does not save at all
    sol=QSS_simGenerate(settings)
    display(plot!(sol[1],sol[2][1]))
    display(plot!(sol[1],sol[2][2]))
    #
end

qssApproachInitialInside()
# @btime qssApproachInitialInside()

 odeDiffEquPackage()

println("done") 
readline()
 =#











#A=[3 -2; 2 -2]; initcond=[1,1]
#anal sol1=(2/3) *exp(2x)+(1/3) *exp(-x)
#analy sol2=(1/3) *exp(2x)+(2/3) *exp(-x)
#= f(x) = (1/3) *exp(2x)+(2/3) *exp(-x)
display(plot!(f, 0, 3))=#
