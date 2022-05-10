 using BenchmarkTools
using StaticArrays
using qss
#= using OrdinaryDiffEq
using DifferentialEquations =#
function qssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0]
    inputVars=@SVector[0.0,1.0]
    jacobian=@SMatrix[0.0 1.0;-1.0 -1.0]
    psettings = ProblemSettings(5.0,saveat(0.1),qss2())
    prob = QSS_Problem(initConditions,jacobian,inputVars)
 
    

    sol=QSS_Solve(psettings,prob)
end

function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2] 
        du[2] = -u[1] - u[2]  
    
    end
    tspan = (0.0,0.10)
    u0 = [1.0,2.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,abstol = 1e-6, reltol = 1e-3)

end

@btime qssApproachInitialInside() 

#@btime odeDiffEquPackage()









#A=[3 -2; 2 -2]; initcond=[1,1]
#anal sol1=(2/3) *exp(2x)+(1/3) *exp(-x)
#analy sol2=(1/3) *exp(2x)+(2/3) *exp(-x)
#= f(x) = (1/3) *exp(2x)+(2/3) *exp(-x)
display(plot!(f, 0, 3))=#
