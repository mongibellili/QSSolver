using BenchmarkTools
using StaticArrays
using qss
#using OrdinaryDiffEq
using DifferentialEquations
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] =        u[2]  - u[3]            
        du[2] = -u[1] -u[2]                       +u[5] 
        du[3] = u[1]         - 2*u[3]      +u[4]      
        du[4] =       u[2]     +u[3]     -2* u[4]  
        du[5] =        2*u[2]   -u[3]               
    end
    tspan = (0.0,10)
    u0 = [1.0,2.0,1.0,-1.0,1.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,abstol = 1e-6, reltol = 1e-3)
    #display(sol)

end
function qssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0,1.0,-1.0,1.0]
    inputVars=@SVector[0.0,0.0,0.0,0.0,0.0]
    jacobian=@SMatrix[0.0 1.0 -1.0 0.0 0.0;
                     -1.0 -1.0 0.0 0.0 1.0;
                      1.0 0.0 -2.0 1.0 0.0;
                      0.0 1.0 1.0 -2.0 0.0;
                      0.0 2.0 -1.0 0.0 0.0]
    #states=2
    psettings = ProblemSettings(5.0,saveat(0.05),qss2())
    prob = QSS_Problem(initConditions,jacobian,inputVars)
 
    

    sol=QSS_Solve(psettings,prob)

    #
end
@btime qssApproachInitialInside()
@btime odeDiffEquPackage()



