#using ModelingToolkit
#using DifferentialEquations
using OrdinaryDiffEq
#using Plots;gr()
using BenchmarkTools
using StaticArrays
#using qss

#using Profile

#using TimerOutputs
#using StaticArrays
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
    tspan = (0.0,1.5)
    u0 = [1.0,2.0]
    prob = ODEProblem(funcName,u0,tspan)
    #,BS3() 
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
    #display(sol)
    display(plot!(sol,line=(:dot, 4)))#, xlims = (0.12,0.15),ylims = (1.23,1.26)))
    #, xlims = (0,0.00005),ylims = (1,1.0001)))
end
function qssApproach()
    initialTime=0.0
    finalTime=0.5 # right now it only accepts float
    dQmin=1e-6
    dQrel=1e-3
    order=1 # order2 means we will consider second derivatives
    solver=qss1()
   # initConditions=[1.0,-1.0,1.0]
   # jacobian=[0.0 1.0 0.0;-1.0 -1.0 0.0;1.0 0.0 -1.0]
    initConditions=[1.0,2.0]
   # jacobianGiven=[0.0 1.0;-1.0 -1.0 ]
    #jacobian=SMatrix{2,2}(jacobianGiven)
    jacobian=@SMatrix[0.0 1.0;-1.0 -1.0 ]
    settings = ModelSettings(initConditions,jacobian,finalTime)
    
   display(settings)
   # simulator = QSS_simulator(settings)
   #QSS_integrate(settings)
    #display(simulator)
    #plotX(simulator)
end
#reset_timer!()
#odeDiffEquPackage()
# qssApproach() 
# @timeit "odeDiffPackage" odeDiffEquPackage()
 #@timeit "qssAproach" qssApproach()
 qssApproach()

#=@profile qssApproach()
f = open("/home/unknown/QS_Solver/prof2.txt", "w")
Profile.print(f)
close(f)
=#
#display(@benchmark qssApproach())

#@time odeDiffEquPackage()
#@btime odeDiffEquPackage()
#readline()
#print_timer()