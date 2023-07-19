using OrdinaryDiffEq
using BSON
using BenchmarkTools
using XLSX
using Plots

function getAverageErrorByRefs(sol::Vector{Vector{Float64}},solRef::Vector{Any},T::Int,numPoints::Int)
  #numPoints=length(solmliqss.savedTimes[1])
  allErrors=0.0
  for index=1:T
      sumTrueSqr=0.0
      sumDiffSqr=0.0
      relerror=0.0
      for i = 1:numPoints #
          ts=solRef[i][index]
          Ns=sol[i][index]
          sumDiffSqr+=(Ns-ts)*(Ns-ts)
          sumTrueSqr+=ts*ts
      end
      relerror=sqrt(sumDiffSqr/sumTrueSqr)
      
      allErrors+= relerror
  end
  return allErrors/T
end

function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[4]-1e6*u[1]+1e3*u[2]
        du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
        du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
        du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
        du[5] = 0.015-200.0*u[2]*u[5]
        du[6] =u[4]-0.6*u[6]
    end
    tspan = (0.0,100.0)
    u0= [0.0,0.75,0.25,0.0,0.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)


 #=    solRodas5P = solve(prob,Rodas5P(),saveat=0.01,abstol = 1e-12, reltol = 1e-8)

    solRodas5PVectorTyson=solRodas5P.u

   
 
   BSON.@save "formalqss/ref_bson/solVectAdvection_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson =#


   absTol=1e-6
   relTol=1e-3

#=  BSON.@load "formalqss/ref_bson/solVectAdvection_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson
 solFTRBDF2 = solve(prob,TRBDF2(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
 err1=getAverageErrorByRefs(solFTRBDF2.u,solRodas5PVectorTyson,6,10000)
 timeTRBDF2=@belapsed solve($prob,TRBDF2(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
 resTRBDF2= ("TRBDF2",err1,timeTRBDF2)
 @show resTRBDF2 =#

 solRosenbrock23 = solve(prob,Rosenbrock23(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
 p1=plot!(solRosenbrock23,marker=(:circle),markersize=2)
 savefig(p1, "plot_solRosenbrock23_TYSON.png")
 
 #= err2=getAverageErrorByRefs(solRosenbrock23.u,solRodas5PVectorTyson,6,10000)
 timeRosenbrock23=@belapsed solve($prob,Rosenbrock23(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
 resRosenbrock23= ("Rosenbrock23",err2,timeRosenbrock23)
 @show resRosenbrock23

 solQBDF2 = solve(prob,QBDF2(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
 err3=getAverageErrorByRefs(solQBDF2.u,solRodas5PVectorTyson,6,10000)
 timeQBDF2=@belapsed solve($prob,QBDF2(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
 resQBDF2= ("QBDF2",err3,timeQBDF2)
 @show resQBDF2
 
 


 XLSX.openxlsx("Tyson diffEqJL solvers_$(relTol)_$(absTol).xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "Tyson diffEqJL solvers__$(relTol)_$(absTol)"
    sheet["A4"] = collect(("solver","error","time"))
    sheet["A5"] = collect(resTRBDF2)
   # sheet["A6"] = collect(resmliqss)
    sheet["A7"] = collect(resRosenbrock23)
    sheet["A8"] = collect(resQBDF2)
 end  =#
end
#@btime 
odeDiffEquPackage()  



