
using qss
using XLSX
#using BenchmarkTools
using BSON

 function tyson(solRodas5PVectorTyson,note,slvr,cycleDetection)
  
    
       odeprob = @NLodeProblem begin
           name=(tyson,)
           u = [0.0,0.75,0.25,0.0,0.0,0.0]
           du[1] = u[4]-1e6*u[1]+1e3*u[2]
           du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
           du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
           du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
           du[5] = 0.015-200.0*u[2]*u[5]
           du[6] =u[4]-0.6*u[6]
       end  
     
    
       println("start solving")
       timenmliqss=0.0
     
       absTol=1e-5
       relTol=1e-2
     
      sol=QSS_Solve(odeprob,slvr,detection=Val(cycleDetection),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=25.0,maxErr=10*relTol)
      solInterp=solInterpolated(sol,0.01)
      err3=getAverageErrorByRefs(solRodas5PVectorTyson,solInterp)
     # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0,maxErr=10*$relTol)
      resnmliqss1E_2= ("$(sol.algName)",relTol,err3,sol.totalSteps,sol.simulStepCount,timenmliqss)
      @show resnmliqss1E_2

      return resnmliqss1E_2
    
end 

BSON.@load "qss/ref_bson/solVect_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson
note="iters_abs-5"
#resmliqss10=tyson(solRodas5PVectorTyson,note,mliqss1_SimulIter(),0)
resmliqss11=tyson(solRodas5PVectorTyson,note,mliqss1_SimulIter(),1)
#= resmliqss12=tyson(solRodas5PVectorTyson,note,mliqss1_SimulIter(),2)
resmliqss13=tyson(solRodas5PVectorTyson,note,mliqss1_SimulIter(),3)
resmliqss14=tyson(solRodas5PVectorTyson,note,mliqss1_SimulIter(),4)

resmliqss20=tyson(solRodas5PVectorTyson,note,mliqss2_SimulIter(),0)
resmliqss21=tyson(solRodas5PVectorTyson,note,mliqss2_SimulIter(),1)
resmliqss22=tyson(solRodas5PVectorTyson,note,mliqss2_SimulIter(),2)
resmliqss23=tyson(solRodas5PVectorTyson,note,mliqss2_SimulIter(),3)
resmliqss24=tyson(solRodas5PVectorTyson,note,mliqss2_SimulIter(),4)

resmliqss30=tyson(solRodas5PVectorTyson,note,mliqss3_SimulIter(),0)
resmliqss31=tyson(solRodas5PVectorTyson,note,mliqss3_SimulIter(),1)
resmliqss32=tyson(solRodas5PVectorTyson,note,mliqss3_SimulIter(),2)
resmliqss33=tyson(solRodas5PVectorTyson,note,mliqss3_SimulIter(),3)
resmliqss34=tyson(solRodas5PVectorTyson,note,mliqss3_SimulIter(),4) =#

#= XLSX.openxlsx("tyson_solvers_$note.xlsx", mode="w") do xf
  sheet = xf[1]
  sheet["A1"] = "tyson_solvers__$note)"
  sheet["A4"] = collect(("solver","Tolerance","error","totalSteps","simul_steps","time"))
  sheet["A5"] = collect(resmliqss10)
  sheet["A6"] = collect(resmliqss11)
  sheet["A7"] = collect(resmliqss12)
  sheet["A8"] = collect(resmliqss13)
  sheet["A9"] = collect(resmliqss14)

  sheet["A10"] = collect(resmliqss20)
  sheet["A11"] = collect(resmliqss21)
  sheet["A12"] = collect(resmliqss22)
  sheet["A13"] = collect(resmliqss23)
  sheet["A14"] = collect(resmliqss24)

  sheet["A15"] = collect(resmliqss30)
  sheet["A16"] = collect(resmliqss31)
  sheet["A17"] = collect(resmliqss32)
  sheet["A18"] = collect(resmliqss33)
  sheet["A19"] = collect(resmliqss34)
 
end  =#