#using ModelingToolkit
using DifferentialEquations
using Plots;gr()
using qss
using BenchmarkTools
using Profile
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
    
        du[1] = u[2]
        du[2] = -u[1] - u[2]
       # du[3] = u[1] - u[3]
    end
    #=@variables x(t) y(t)
    @derivatives D'~t
    eqs = [D(x) ~ (x),D(y) ~ -x-y]
    de = ODESystem(eqs)
    funcName = ODEFunction(de, [x,y])
    =#
    tspan = (0.0,5)
    u0 = [1.0,2.0]
    prob = ODEProblem(funcName,u0,tspan)
    #,BS3() 
    sol = solve(prob,abstol = 1e-6, reltol = 1e-3)
    #display(sol)
    display(plot!(sol,line=(:dot, 4)))#, xlims = (0.12,0.15),ylims = (1.23,1.26)))
    #, xlims = (0,0.00005),ylims = (1,1.0001)))
end
function qssApproach()
    initialTime=0.0
    finalTime=2.5 # right now it only accepts float
    dQmin=1e-6
    dQrel=1e-3
    order=1 # order2 means we will consider second derivatives
    solver="qss1"
   # initConditions=[1.0,-1.0,1.0]
   # jacobian=[0.0 1.0 0.0;-1.0 -1.0 0.0;1.0 0.0 -1.0]
    initConditions=[1.0,2.0]
    jacobian=[0.0 1.0;-1.0 -1.0 ]

    settings = ModelSettings(initialTime,finalTime,dQmin,dQrel,order,solver,initConditions,jacobian)
    

    simulator = QSS_simulator(settings)
    QSS_simulate(simulator)
    #display(simulator.qssData.x)
    #plotX(simulator)
end
#odeDiffEquPackage()
#qssApproach()
@btime qssApproach()
#=@profile qssApproach()
f = open("/home/unknown/QS_Solver/prof2.txt", "w")
Profile.print(f)
close(f)
=#
#display(@benchmark qssApproach())

#@btime odeDiffEquPackage()
#display(@benchmark odeDiffEquPackage())
#readline()