using qss

#include("D://QS_Solver//qss//src//models//classicProblem.jl") 
include("/home/unknown/QS_Solver/qss/src/models/classicProblem.jl")# where you saved the model
function test()
    odeprob = @NLodeProblem begin
        
        #----------twoVarSys1-------------
        
       #=  u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]  =#
       
        #------------twoVarSys12-----------
        #= du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0 =#

        #------------twoVarSys13----------
        u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2
        #------------twoVarSys14----------
        #= du[1] = -u[1]-10.0*u[2]+0.2
        du[2] =10.0*u[1]-u[2]+1.2 =#
      
    end
     
   # sol=QSS_Solve_from_model(twoVarSys13,odeprob,15.0,liqss2(),saveat(0.00001)) # 
  #  sol=QSS_Solve(odeprob,10.0,liqss2())
    sol1=QSS_Solve(odeprob,5.0,mliqss2(),saveat(0.01)) 
   # sol2=QSS_Solve(odeprob,50.0,mliqss1(),saveat(0.001)) 

   # sol=QSS_Solve(odeprob,5.0,liqss2())
  # display(sol(1,14.0))
 #plotSol(sol1)
 #plotSol(sol2)

end
test()

