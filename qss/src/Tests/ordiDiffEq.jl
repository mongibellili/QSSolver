#using OrdinaryDiffEq
using DifferentialEquations
#using qssv01
using BenchmarkTools
using Plots
function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
        #= du[1] = -5*u[1]
        du[2] =-u[1]-5*u[2] =#
        #= du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2 =#
        du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0

        #= du[1] = -4.0*u[1]+0.5*u[2]+1.0
        du[2] =u[1]-0.25*u[2]+0.25 =#
      #=   du[1] = u[2]
        du[2] =-u[1]-u[2]  =#
    
    end
    tspan = (0.0,1000.0)
    u0 = [1.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),saveat=1.0,abstol = 1e-9, reltol = 1e-5)
    #sol = solve(prob,Rosenbrock23(),abstol = 1e-6, reltol = 1e-3)
  
    u1=(7*sqrt(51)-50)/100
    u2=(-7*sqrt(51)-50)/100
    λ1=-7*sqrt(51)-50
    λ2=7*sqrt(51)-50
    c2=960/(7*sqrt(51))
    c1=-c2
    
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+20.2
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)
   
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
    display(relerror1)  =#

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
odeDiffEquPackage()  

#=u1=(-sqrt(257)-15)/8
u2=(sqrt(257)-15)/8
λ1=(-sqrt(257)-17)/8
λ2=(sqrt(257)-17)/8
c1=-4/sqrt(257)
c2=4/sqrt(257)
x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)
x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)
display(plot!(x1,title="against BS3:ΔQ=1e-5",label="true x1",xlims=(100,160),ylims=(-0.000002,0.000002)))
display(plot!(x2,label="true x2",xlims=(100,160),ylims=(-0.000002,0.000002)))
println("press enter to exit")
readline()
 =#

#= function test()
    odeprob = @NLodeProblem begin
        u = [1.0, 0.0]
        discrete = [0.0]
        #=         du[1] = u[2]
        du[2] =-u[1]-u[2] =#
        du[1] = 0.01*u[2]
        du[2] =-100.0*u[1]-100.0*u[2]+2020.0
        #= du[1] = -u[1]-u[2]+0.2
        du[2] =u[1]-u[2]+1.2 =#
    end
   sol = QSS_Solve(odeprob,2.0,liqss2())
   plotSol(sol)
 
   
end
test() =#
#odeDiffEquPackage() 