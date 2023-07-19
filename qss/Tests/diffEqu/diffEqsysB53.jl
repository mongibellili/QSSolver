using OrdinaryDiffEq
using BenchmarkTools
using XLSX
using Plots

function getAverageErrorByRefs(sol::Vector{Vector{Float64}},solT::Vector{Float64},solRef::Function,numPoints::Int,index::Int)
  #numPoints=length(solmliqss.savedTimes[1])
  
      sumTrueSqr=0.0
      sumDiffSqr=0.0
      relerror=0.0
      for i = 1:numPoints #
          ts=solRef(solT[i])
          Ns=sol[i][index]
          sumDiffSqr+=(Ns-ts)*(Ns-ts)
          sumTrueSqr+=ts*ts
      end
      relerror=sqrt(sumDiffSqr/sumTrueSqr)
      
     
  
  return relerror
end


function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
        du[1] = -20.0*u[1]-80.0*u[2]+1600.0
         du[2] =1.24*u[1]-0.01*u[2]+0.2
    end
    tspan = (0.0,10.0)
    u0= [-1.0, -2.0]
    prob = ODEProblem(funcName,u0,tspan)


    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
    c1, c2 = 121.14809142478035, -143.14809142478035
    xp1, xp2 = 0.0, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 


   absTol=1e-5
   relTol=1e-2

#= 
 solFTRBDF2 = solve(prob,TRBDF2(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
 err1=getAverageErrorByRefs(solFTRBDF2.u,solFTRBDF2.t,x1,10000,1)
 err2=getAverageErrorByRefs(solFTRBDF2.u,solFTRBDF2.t,x2,10000,2)
 timeTRBDF2=@belapsed solve($prob,TRBDF2(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
 resTRBDF2= ("TRBDF2",(err1+err2)/2,timeTRBDF2)
 @show resTRBDF2 =#

 solRosenbrock23 = solve(prob,Rosenbrock23(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
 p1=plot!(solRosenbrock23,marker=(:circle),markersize=2)
 savefig(p1, "plot_solRosenbrock23_sys53.png")
 #= err1=getAverageErrorByRefs(solRosenbrock23.u,solRosenbrock23.t,x1,10000,1)
 err2=getAverageErrorByRefs(solRosenbrock23.u,solRosenbrock23.t,x2,10000,2)
 timeRosenbrock23=@belapsed solve($prob,Rosenbrock23(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
 resRosenbrock23= ("Rosenbrock23",(err1+err2)/2,timeRosenbrock23)
 @show resRosenbrock23 =#

 #= solQBDF2 = solve(prob,QBDF2(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
 err1=getAverageErrorByRefs(solQBDF2.u,solQBDF2.t,x1,10000,1)
 err2=getAverageErrorByRefs(solQBDF2.u,solQBDF2.t,x2,10000,2)
 timeQBDF2=@belapsed solve($prob,QBDF2(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
 resQBDF2= ("QBDF2",(err1+err2)/2,timeQBDF2)
 @show resQBDF2
 
 


 XLSX.openxlsx("sysB53 diffEqJL solvers_$(relTol)_$(absTol).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "sysB53 diffEqJL solvers__$(relTol)_$(absTol)"
    sheet["A4"] = collect(("solver","error","time"))
    sheet["A5"] = collect(resTRBDF2)
   # sheet["A6"] = collect(resmliqss)
    sheet["A7"] = collect(resRosenbrock23)
    sheet["A8"] = collect(resQBDF2)
 end  =#
end
#@btime 
odeDiffEquPackage()  



