
using qss
using XLSX
#using BenchmarkTools
using BSON

 function Oregon(solRodas5PVectorOregon,note,slvr,cycleDetection)
  
    
       odeprob = @NLodeProblem begin
        name=(oregonator,)
    u  = [1.0,1.0,0.0]
    du[1] = 100.8*(9.523809523809524e-5*u[2]-u[1]*u[2]+u[1]*(1.0-u[1]))
    du[2] =40320.0*(-9.523809523809524e-5*u[2]-u[1]*u[2]+u[3])
    du[3] = u[1] -u[3]
       end  
     
    
       println("start solving")
       timenmliqss=0.0
     
       absTol=1e-3
       relTol=1e-2
     
      sol=QSS_Solve(odeprob,slvr,detection=Val(cycleDetection),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=25.0#= ,maxErr=10*relTol =#)
     # save_Sol(sol)
      solInterp=solInterpolated(sol,0.01)
      err3=0.0
      err3=getAverageErrorByRefs(solRodas5PVectorOregon,solInterp)
     # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0,maxErr=10*$relTol)
      resnmliqss1E_2= ("$(sol.algName)",err3,sol.totalSteps,sol.simulStepCount,timenmliqss)
      @show resnmliqss1E_2

      return resnmliqss1E_2
    
end 

BSON.@load "qss/ref_bson/solVectAdvection_Oregon_Rodas5Pe-12.bson" solRodas5PVectorOregon
#BSON.@load "ref_bson/solVect_Tyson_Rodas5Pe-12.bson" solRodas5PVectorOregon
note="Ana_abs-3_1Delta"
resmliqss10=Oregon(solRodas5PVectorOregon,note,mliqss1_SimulAna(),0)
resmliqss11=Oregon(solRodas5PVectorOregon,note,mliqss1_SimulAna(),1)
resmliqss12=Oregon(solRodas5PVectorOregon,note,mliqss1_SimulAna(),2)
resmliqss13=Oregon(solRodas5PVectorOregon,note,mliqss1_SimulAna(),3)
resmliqss14=Oregon(solRodas5PVectorOregon,note,mliqss1_SimulAna(),4)

resmliqss20=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulAna(),0)
resmliqss21=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulAna(),1)
resmliqss22=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulAna(),2)
resmliqss23=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulAna(),3)
resmliqss24=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulAna(),4)

resmliqss30=Oregon(solRodas5PVectorOregon,note,mliqss3_SimulAna(),0)
resmliqss31=Oregon(solRodas5PVectorOregon,note,mliqss3_SimulAna(),1)
resmliqss32=Oregon(solRodas5PVectorOregon,note,mliqss3_SimulAna(),2)
resmliqss33=Oregon(solRodas5PVectorOregon,note,mliqss3_SimulAna(),3)
resmliqss34=Oregon(solRodas5PVectorOregon,note,mliqss3_SimulAna(),4)

XLSX.openxlsx("Oregon_solvers_$note.xlsx", mode="w") do xf
  sheet = xf[1]
  sheet["A1"] = "Oregon_solvers__$note)"
  sheet["A4"] = collect(("solver","error","totalSteps","simul_steps","time"))
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
 
end 