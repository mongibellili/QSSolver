using qss
using BenchmarkTools
using Plots
#include("D://QS_Solver//qss//src//models//classicProblem.jl") 
include("/home/unknown/QS_Solver/qss/src/models/classicProblem.jl")# where you saved the model
function test()
    odeprob = @NLodeProblem begin
        
        #----------twoVarSys1-------------
        
        #= u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]  =#
       
        #------------twoVarSys12-----------
        
        #= u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0 =#

        #------------twoVarSys13----------
        #= u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2 =#
        #------------twoVarSys14----------
        #= u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-10.0*u[2]+0.2
        du[2] =10.0*u[1]-u[2]+1.2 =#
        #------------twoVarSys15----------
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = -4.0*u[1]+0.5*u[2]
        du[2] =u[1]-0.25*u[2] 
        #------------twoVarSys16----------
        #= u = [1.0, 2.3]
        discrete = [0.0]
        du[1] = -1.01*u[1]-100.0*u[2]+1000.2
        du[2] =10.1*u[1]-1.001*u[2]+0.2 =#
      
    end
    

    #= solqss1=QSS_Solve_from_model(twoVarSys15,odeprob,100.0,qss1(),saveat(0.01),0.0,1e-6,1e-3)    
    solqss2=QSS_Solve_from_model(twoVarSys15,odeprob,100.0,qss2(),saveat(0.01),0.0,1e-6,1e-3) # 
    solqss3=QSS_Solve_from_model(twoVarSys15,odeprob,100.0,qss3(),saveat(0.01),0.0,1e-6,1e-3)
    solliqss1=QSS_Solve_from_model(twoVarSys15,odeprob,100.0,liqss1(),saveat(0.01),0.0,1e-6,1e-3)
    solliqss2=QSS_Solve_from_model(twoVarSys15,odeprob,100.0,liqss2(),saveat(0.01),0.0,1e-6,1e-3)
    solliqss3=QSS_Solve_from_model(twoVarSys15,odeprob,100.0,liqss3(),saveat(0.01),0.0,1e-6,1e-3)
    solmliqss1=QSS_Solve_from_model(twoVarSys15,odeprob,100.0,mliqss1(),saveat(0.01),0.0,1e-6,1e-3) =#
    sol=QSS_Solve_from_model(twoVarSys15,odeprob,50.0,qss1(),saveat(0.1),0.0,1e-9,1e-5)
    

     u1=(-sqrt(257)-15)/8
    u2=(sqrt(257)-15)/8
    λ1=(-sqrt(257)-17)/8
    λ2=(sqrt(257)-17)/8
    c1=-4/sqrt(257)
    c2=4/sqrt(257)
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)
    plotError(sol,1,x1)
#=
    display(getError(solqss1,1,x1));println()
    display(getError(solqss2,1,x1));println()
    display(getError(solqss3,1,x1));println()
    display(getError(solliqss1,1,x1));println()
    display(getError(solliqss2,1,x1));println()
    display(getError(solliqss3,1,x1));println()
    display(getError(solmliqss1,1,x1));println()
    display(getError(solmliqss2,1,x1));println()
    println(" variable 2: ")
    display(getError(solqss1,2,x2));println()
    display(getError(solqss2,2,x2));println()
    display(getError(solqss3,2,x2));println()
    display(getError(solliqss1,2,x2));println()
    display(getError(solliqss2,2,x2));println()
    display(getError(solliqss3,2,x2));println()
    display(getError(solmliqss1,2,x2));println()
    display(getError(solmliqss2,2,x2));println()
     =#



   #sol=QSS_Solve(odeprob,20.0,mliqss2())
  # display(sol(1,14.0))
 #plotSol(sol)
 #plotSol(sol2)
#=  display(evaluateSol(sol,1,0.3));println()
 display(evaluateSol(sol,1,0.4));println()
 display(evaluateSol(sol,1,0.5));println()
 display(evaluateSol(sol,2,0.2));println()
 display(evaluateSol(sol,2,0.3));println()
 display(evaluateSol(sol,2,0.4));println()
 display(evaluateSol(sol,2,0.5));println() =#

end
#@btime 
test()

#= 
u1=(-sqrt(257)-15)/8
u2=(sqrt(257)-15)/8
λ1=(-sqrt(257)-17)/8
λ2=(sqrt(257)-17)/8
c1=-4/sqrt(257)
c2=4/sqrt(257)
x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)
x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)
display(plot!(x1,title="against liqss3: quan-5  ",label="true x1",xlims=(100,160),ylims=(-0.000002,0.000002)))
display(plot!(x2,label="true x2",xlims=(100,160),ylims=(-0.000002,0.000002)))
println("press enter to exit")
readline() =#