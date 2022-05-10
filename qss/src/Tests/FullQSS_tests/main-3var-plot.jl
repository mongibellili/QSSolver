
using Plots;gr()

using StaticArrays
using qss
using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2] - u[3] 
        du[2] = -u[1] - u[2]  
        du[3] = u[1] - 2*u[3] 
    end
    tspan = (0.0,5)
    u0 = [1.0,2.0,1.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end
function qssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0,1.0]
    jacobian=@SMatrix[0.0 1.0 -1.0;-1.0 -1.0 0.0;1.0 0.0 -2.0]
    #states=2
    psettings = ProblemSettings(5.0,saveat(0.1),qss1())
    prob = QSS_Problem(initConditions,jacobian,inputVars)
 
    

    sol=QSS_Solve(psettings,prob)
    display(plot!(sol[1],sol[2][1]))
    display(plot!(sol[1],sol[2][2]))
    display(plot!(sol[1],sol[2][3]))
    #
end
qssApproachInitialInside()
odeDiffEquPackage()
println("done") 
readline()


