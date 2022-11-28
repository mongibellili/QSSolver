using qss
using BenchmarkTools
#using Plots
#using OrdinaryDiffEq

#include("D://QS_Solver//qss//src//models//classicProblem.jl") 
include("/home/unknown/QS_Solver/qss/src/models/classicProblem.jl")# where you saved the model
function test()
    odeprob = @NLodeProblem begin
        
        #----------twoVarSys1-------------
        
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2] 
       
        #------------twoVarSys12-----------
        
       #=  u = [1.0, 0.0]
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
    

    




    #solqss1=QSS_Solve_from_model(twoVarSys15,odeprob,50.0,qss1(),saveat(0.5),0.0,1e-9,1e-5)    
    solqss2=QSS_Solve_from_model(twoVarSys1,odeprob,5.0,qss2(),saveat(0.1),0.0,1e-6,1e-3) 
    solqss3=QSS_Solve_from_model(twoVarSys1,odeprob,5.0,qss3(),saveat(0.1),0.0,1e-6,1e-3)
   
    #solliqss1=QSS_Solve_from_model(twoVarSys15,odeprob,50.0,liqss1(),saveat(0.5),0.0,1e-9,1e-5)
    solliqss2=QSS_Solve_from_model(twoVarSys1,odeprob,5.0,liqss2(),saveat(0.1),0.0,1e-6,1e-3)
    solliqss3=QSS_Solve_from_model(twoVarSys1,odeprob,5.0,liqss3(),saveat(0.1),0.0,1e-6,1e-3)

    #solmliqss1=QSS_Solve_from_model(twoVarSys15,odeprob,50.0,mliqss1(),saveat(0.5),0.0,1e-9,1e-5)
    solmliqss2=QSS_Solve_from_model(twoVarSys1,odeprob,5.0,mliqss2(),saveat(0.1),0.0,1e-6,1e-3)
    

   #=  u1=(-sqrt(257)-15)/8
    u2=(sqrt(257)-15)/8
    λ1=(-sqrt(257)-17)/8
    λ2=(sqrt(257)-17)/8
    c2=(397-2*(sqrt(257)+15))/sqrt(257)
    c1=-4-c2
    b1=1.0
    b2=0.25
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+0.5*b1+b2
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+2*b1+8*b2 =#
#=     plotError(solqss1,1,x1)
    plotError(solqss2,1,x1)
    plotError(solqss3,1,x1)
    plotError(solliqss1,1,x1)
    plotError(solliqss2,1,x1)
    plotError(solliqss3,1,x1)
    plotError(solmliqss1,1,x1)
    plotError(solmliqss2,1,x1) =#
  #  plotError(solqss1,2,x2)
   # plotError(solqss2,2,x2)
   # plotError(solqss3,2,x2)
   # plotError(solliqss1,2,x2)
  #  plotError(solliqss2,2,x2)
  #  plotError(solliqss3,2,x2)
   # plotError(solmliqss1,2,x2)
   # plotError(solmliqss2,2,x2)
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



  # sol=QSS_Solve(odeprob,20.0,liqss2())
  # sol=QSS_Solve(odeprob,50.0,qss3(),saveat(0.1),0.0,1e-6,1e-3)
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
function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
        #= du[1] = -5*u[1]
        du[2] =-u[1]-5*u[2] =#
        #= du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2 =#
        #= du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0 =#

        du[1] = -4.0*u[1]+0.5*u[2]+1.0
        du[2] =u[1]-0.25*u[2]+0.25
    
    end
    tspan = (0.0,50.0)
    u0 = [100.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),saveat=0.5,abstol = 1e-9, reltol = 1e-5)
    #sol = solve(prob,Rosenbrock23(),abstol = 1e-6, reltol = 1e-3)
  
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
   #= sumTrueSqr1=0.0
   sumDiffSqr1=0.0
   relerror1=0.0
   numPoints=length(sol.t)
   for i = 1:numPoints
      ft=x2(sol.t[i])
      sumDiffSqr1+=(sol.u[i][2]-ft)*(sol.u[i][2]-ft)
      sumTrueSqr1+=ft*ft
    end
    relerror1=sqrt(sumDiffSqr1/sumTrueSqr1)
    display(relerror1) =#

    numPoints=length(sol.t)
   
    
      temp = []
      for i = 1:numPoints #each point is a taylor
        ft=x2(sol.t[i])
        push!(temp, abs(sol.u[i][2]-ft)/ft)
      end
     display(plot!(sol.t, temp,label="BS3")) 
     println("done")
     readline()

   # display(plot!(sol))
 #=    display(plot!(sol,xlims=(100,160),ylims=(-0.000002,0.000002)))
  #  display(plot!(sol,xlims=(10.0,10.0005),ylims=(-0.49971,-0.49960)))
    println("done")
    readline()  =#
end
#@btime 

odeDiffEquPackage()   =#









#= 
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
#= display(plot!(x1,title="against liqss3: quan-5  ",label="true x1",xlims=(0,1),ylims=(-0.000002,12.000002)))
display(plot!(x2,label="true x2",xlims=(0,1),ylims=(-0.000002,12.000002))) =#
display(plot!(x1,title="against liqss3: quan-5  ",label="true x1",xlims=(0,20)))
display(plot!(x2,label="true x2",xlims=(0,20)))
println("press enter to exit")
readline() =#