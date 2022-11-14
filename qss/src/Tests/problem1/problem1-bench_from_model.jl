using qss
using BenchmarkTools
include("/home/unknown/QS_Solver/qss/src/models/classicProblem.jl") 


function test()
     odeprob = @NLodeProblem begin
       
      #------------twoVarSys1----------- 
      #= u = [1.0, 0.0]
          discrete = [0.0]
           du[1] = u[2]
          du[2] =-u[1]-u[2]  =#

           #------------twoVarSys12-----------
           #= u = [1.0, 0.0]
          discrete = [0.0]
        du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0 =#
        #------------twoVarSys13-----------
        u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2

      end
     sol=QSS_Solve_from_model(twoVarSys13,odeprob,30.0,liqss3())
end
@btime test()
