using Plots;gr()

using StaticArrays
using qss
using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2] #+t+1
        du[2] = -u[1] - u[2]  +1#+10-t*t+t +cos(t)
    
    end
    tspan = (0.0,5)
    u0 = [0.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end

function qssApproachInitialInside() 
    initConditions=@SVector[0.0,0.0]
    jacobian=@SMatrix[0.0 1.0;-1.0 -1.0 ]
 
#=     u1(t)=t+1   
    u2(t)=10-t*t+t #+cos(t)
    inputVars=SVector{2,Function}(u1,u2) =#   # later a multiplexer to map each inputvar to its correspondant state var

    inputVars=@SVector[0.0,1.0]

    psettings = ProblemSettings(5.0,saveat(0.1),qss2())
    prob = QSS_Problem(initConditions,jacobian,inputVars)
 
    

    sol=QSS_Solve(psettings,prob)
    display(plot!(sol[1],sol[2][1]))
    display(plot!(sol[1],sol[2][2])) 
   
end

qssApproachInitialInside()


 odeDiffEquPackage()

println("done") 
readline()








