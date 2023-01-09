using qss
using BenchmarkTools
using Plots
#using OrdinaryDiffEq
#using DifferentialEquations

include("D://QS_Solver//qss//src//models//classicProblem.jl") 
#include("/home/unknown/QS_Solver/qss/src/models/classicProblem.jl")# where you saved the model
function test()
    odeprob = @NLodeProblem begin 
       
      
       #------------twoVarSys12-----------
       
       u = [1.0, 0.0]
       discrete = [0.0]
       du[1] = 0.01*u[2]
       du[2] =-100.0*u[1]-100.0*u[2]+2020.0

     
    end
   # solliqss1=QSS_Solve_from_model(twoVarSys12,odeprob,50.0,liqss1(),saveat(0.5),0.0,1e-9,1e-5)
<<<<<<< HEAD
    solliqss2=QSS_Solve_from_model(twoVarSys12,odeprob,1000.0,liqss2(),saveat(1.0),0.0,1e-9,1e-5)
    solliqss3=QSS_Solve_from_model(twoVarSys12,odeprob,1000.0,liqss3(),saveat(1.0),0.0,1e-9,1e-5)
   # solliqss3=QSS_Solve_from_model(twoVarSys12,odeprob,5.0,liqss3())

   u1=(7*sqrt(51)-50)/100
=======
   # solliqss2=QSS_Solve_from_model(twoVarSys12,odeprob,200.0,liqss2(),saveat(0.5),0.0,1e-9,1e-5)
    solliqss3=QSS_Solve_from_model(twoVarSys12,odeprob,5.0,liqss3(),saveat(0.1),0.0,1e-6,1e-3)
   # solliqss3=QSS_Solve_from_model(twoVarSys12,odeprob,5.0,liqss3())

   #= u1=(7*sqrt(51)-50)/100
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
   u2=(-7*sqrt(51)-50)/100
   λ1=-7*sqrt(51)-50
   λ2=7*sqrt(51)-50
   c2=960/(7*sqrt(51))
   c1=-c2
   
   x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+20.2
   x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)
<<<<<<< HEAD
     plotRelativeError(solliqss2,1,x1)
     plotRelativeError(solliqss3,1,x1)
   #  plotAbsoluteError(solliqss3,2,x2)
    #plotSquaredRelativeError(solliqss2,2,x2)
   #=  plotMSE(solliqss2,2,x2) 
    plotMSE(solliqss3,2,x2)  =#

=======
    # plotError(solliqss1,2,x2)
    plotRelativeError(solliqss2,2,x2)
    plotRelativeError(solliqss3,2,x2)  =#
    #= plotCumulativeSquaredRelativeError(solliqss2,2,x2)
    plotCumulativeSquaredRelativeError(solliqss3,2,x2)  =#
    
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3

   # plotSol(solliqss2)
 #plotSol(solliqss3)
end
@btime test()

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
  # sol = solve(prob,BS3()) 
   sol = solve(prob,Rodas5P()) 


 display(plot!(sol))
 println("done")
 readline()
 
end
 odeDiffEquPackage()  =#  
<<<<<<< HEAD

#=  u1=(7*sqrt(51)-50)/100
=======
#= 
 u1=(7*sqrt(51)-50)/100
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
u2=(-7*sqrt(51)-50)/100
λ1=-7*sqrt(51)-50
λ2=7*sqrt(51)-50
c2=960/(7*sqrt(51))
c1=-c2

x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+20.2
x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)
#= display(plot!(x1,title="against liqss3: quan-5  ",label="true x1",xlims=(0,1),ylims=(-0.000002,12.000002)))
display(plot!(x2,label="true x2",xlims=(0,1),ylims=(-0.000002,12.000002))) =#
#display(plot!(x1,title="sys2  true  ",label="true x1",xlims=(0,5)))
display(plot!(x1,label="true x1",xlims=(0,1000)))
display(plot!(x2,label="true x2",xlims=(0,1000)))
println("press enter to exit")
readline() =#