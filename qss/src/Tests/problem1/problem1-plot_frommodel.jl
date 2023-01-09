using qss
using BenchmarkTools
using Plots
#using OrdinaryDiffEq
#using DifferentialEquations

include("D://QS_Solver//qss//src//models//classicProblem.jl") 
#include("/home/unknown/QS_Solver/qss/src/models/classicProblem.jl")# where you saved the model
function test()
    odeprob = @NLodeProblem begin 
       #----------twoVarSys1-------------
        
       #= u = [1.0, 0.0]
         discrete = [0.0]
         du[1] = u[2]
         du[2] =-u[1]-u[2]  =#
      
       #------------twoVarSys12-----------
       
    #=    u = [1.0, 0.0]
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
         u = [100.0, 0.0]
         discrete = [0.0]
         du[1] = -4.0*u[1]+0.5*u[2]+1.0
         du[2] =u[1]-0.25*u[2]+0.25
         #------------twoVarSys16----------
         #= u = [1.0, 2.3]
         discrete = [0.0]
         du[1] = -1.01*u[1]-100.0*u[2]+1000.2
         du[2] =10.1*u[1]-1.001*u[2]+0.2 =#
     
    end
   # solliqss1=QSS_Solve_from_model(twoVarSys15,odeprob,50.0,liqss1(),saveat(0.5),0.0,1e-9,1e-5)
     solliqss2=QSS_Solve_from_model(twoVarSys15,odeprob,1000.0,liqss2(),saveat(1.5),0.0,1e-9,1e-5)
    solliqss3=QSS_Solve_from_model(twoVarSys15,odeprob,1000.0,liqss3(),saveat(1.5),0.0,1e-9,1e-5)
  # solliqss3=QSS_Solve_from_model(twoVarSys12,odeprob,5.0,liqss3())

    u1=(-sqrt(257)-15)/8
    u2=(sqrt(257)-15)/8
    λ1=(-sqrt(257)-17)/8
    λ2=(sqrt(257)-17)/8
    c2=(397-2*(sqrt(257)+15))/sqrt(257)
    c1=-4-c2
    b1=1.0
    b2=0.25
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+0.5*b1+b2
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+2*b1+8*b2 
    # plotError(solliqss1,2,x2)
    plotError(solliqss2,2,x2)
    plotError(solliqss3,2,x2)


 #plotSol(solliqss3)
end
#@btime 
test()

#= function odeDiffEquPackage()
 
   function funcName(du,u,p,t)# api requires four args
       #= du[1] = -5*u[1]
       du[2] =-u[1]-5*u[2] =#
       #= du[1] = -u[1]-u[2]+0.2
       du[2] =u[1]-u[2]+1.2 =#
       #= du[1] = 0.01*u[2]
       du[2] =-100.0*u[1]-100.0*u[2]+2020.0 =#

       du[1] = 0.01*u[2]
       du[2] =-100.0*u[1]-100.0*u[2]+2020.0
   
   end
   tspan = (0.0,5.0)
   u0 = [1.0,0.0]
   prob = ODEProblem(funcName,u0,tspan)
   sol = solve(prob,BS3())


 display(plot!(sol))
 println("done")
 readline()
 
end
   #@btime 
 odeDiffEquPackage()   
 =#