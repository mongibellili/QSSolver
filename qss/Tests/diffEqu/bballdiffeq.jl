#using OrdinaryDiffEq # stuck with events
using DifferentialEquations
#using qssv01
#using BenchmarkTools
using Plots
function odeDiffEquPackage()
    function f(du,u,p,t)
        du[1] = u[2]
        du[2] = -9.8
    end
    function condition(u,t,integrator) # Event when event_f(u,t) == 0
        u[1]
    end
    function affect!(integrator)
        integrator.u[2] = -integrator.u[2]
    end
   
    cb = ContinuousCallback(condition,affect!)
   
    u0 = [10.0,0.0]
    tspan = (0.0,2.0)
    p = 9.8
    prob = ODEProblem(f,u0,tspan)
  
    sol = solve(prob,Tsit5(),callback=cb)
   
    p1=plot!(sol,xlims=(1.41,1.45) ,ylims=(-2.04e-1,2.06e-1))
    savefig(p1, "tsit5_bball")
 end
odeDiffEquPackage() 