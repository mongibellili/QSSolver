using qss
using BenchmarkTools
#using Plots
#using OrdinaryDiffEq

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
       
       u = [1.0, 0.0]
       discrete = [0.0]
       du[1] = 0.01*u[2]
       du[2] =-100.0*u[1]-100.0*u[2]+2020.0

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
       #= u = [100.0, 0.0]
       discrete = [0.0]
       du[1] = -4.0*u[1]+0.5*u[2]+1.0
       du[2] =u[1]-0.25*u[2]+0.25 =#
       #------------twoVarSys16----------
       #= u = [1.0, 2.3]
       discrete = [0.0]
       du[1] = -1.01*u[1]-100.0*u[2]+1000.2
       du[2] =10.1*u[1]-1.001*u[2]+0.2 =#
     
    end
    solmliqss2=QSS_Solve_from_model(twoVarSys12,odeprob,5.0,liqss3())
 #plotSol(solmliqss2)
end
#@btime 
test()
