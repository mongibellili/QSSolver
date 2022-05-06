using Plots;gr()

using StaticArrays
using qss
#= using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2] 
        du[2] = -u[1] - u[2]  
    
    end
    tspan = (0.0,5)
    u0 = [1.0,2.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end =#

function qssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0]
    jacobian=@SMatrix[0.0 1.0;-1.0 -1.0 ]
    #states=2
    settings = ModelSettings(initConditions,jacobian,5.7,saveat(0.1),qss2())#do not call saveat to not save, i should fix when called with zero; it does not save at all
    sol=QSS_simGenerate(settings)
    display(plot!(sol[1],sol[2][1]))
    display(plot!(sol[1],sol[2][2]))
    #
end

qssApproachInitialInside()


 #odeDiffEquPackage()

println("done") 
readline()








