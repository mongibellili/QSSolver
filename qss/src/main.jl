
using Plots;gr()
using BenchmarkTools
using StaticArrays
using qss
using OrdinaryDiffEq
#using TimerOutputs  
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2]
        du[2] = -u[1] - u[2]   
    end
    tspan = (0.0,3)
    u0 = [1.0,2.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
    display(plot!(sol,line=(:dot, 4)))
end
function qssApproach(initConditions::Vector{Float64},jacobian::Matrix{Float64}) 
    settings = ModelSettings(initConditions,jacobian,0.5,saveat(0.15))
    sol=QSS_integrate(settings)
    #display(plot!(sol[1],sol[2][1]))
end
function qssApproach2(initConditions::Vector{Float64},jacobian::Matrix{Float64}) 
    settings = ModelSettings(initConditions,jacobian,3.0,0.0,1e-4,1e-2,qss1(),saveat(0.05))
    sol=QSS_integrate(settings)
    display(plot!(sol[1],sol[2][1],line=(:dash, 2)))
end
initConditions=[1.0,2.0]
jacobian=[0.0 1.0;-1.0 -1.0]
#= @btime qssApproach(initConditions,jacobian) 
@btime qssApproach2(initConditions,jacobian)  =#

qssApproach(initConditions,jacobian) 
#qssApproach2(initConditions,jacobian) 
# @btime qssApproach(initConditions,jacobian)
# odeDiffEquPackage()
#A=[3 -2; 2 -2]; initcond=[1,1]
#anal sol1=(2/3) *exp(2x)+(1/3) *exp(-x)
#analy sol2=(1/3) *exp(2x)+(2/3) *exp(-x)



#= f(x) = (1/3) *exp(2x)+(2/3) *exp(-x)
display(plot!(f, 0, 3))=#
#println("done") 
readline()

