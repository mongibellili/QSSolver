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
      #=   u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2 =#
        #------------twoVarSys14----------
        #= du[1] = -u[1]-10.0*u[2]+0.2
        du[2] =10.0*u[1]-u[2]+1.2 =#
        #------------twoVarSys15----------
        u = [0.5, 0.0]
        discrete = [0.0]
        du[1] = -4.0*u[1]+0.5*u[2]
        du[2] =u[1]-0.25*u[2] 
      
    end
     
    sol1=QSS_Solve_from_model(twoVarSys15,odeprob,160.0,qss3(),saveat(0.001),0.0,1e-9,1e-5) # 
   # sol1=QSS_Solve(odeprob,40.0,liqss2())
    #sol1=QSS_Solve(odeprob,30.0,mliqss1(),saveat(0.01)) 
   # sol2=QSS_Solve_from_model(twoVarSys13,odeprob,50.0,mliqss1(),saveat(0.01)) # 

   # sol=QSS_Solve(odeprob,5.0,liqss2())
  # display(sol(1,14.0))
plotSol(sol1)
 #plotSol(sol2)

end
#@btime 
test()


u1=(-sqrt(257)-15)/8
u2=(sqrt(257)-15)/8
λ1=(-sqrt(257)-17)/8
λ2=(sqrt(257)-17)/8
c1=-4/sqrt(257)
c2=4/sqrt(257)
x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)
x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)
display(plot!(x1,title="against qss3:ΔQ=1e-5",label="true x1",xlims=(100,160),ylims=(-0.000002,0.000002)))
display(plot!(x2,label="true x2",xlims=(100,160),ylims=(-0.000002,0.000002)))
println("press enter to exit")
readline()