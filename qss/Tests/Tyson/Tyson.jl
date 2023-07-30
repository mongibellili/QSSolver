

using qss
using XLSX
using BenchmarkTools
using BSON

function test53()
  BSON.@load "qss/ref_bson/solVectAdvection_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson
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

     absTol=1e-6
     relTol=1e-3

   

   #  solliqss2=QSS_Solve(odeprob,liqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=100.0)
    # save_Sol(solliqss2)
   # solmliqss2=QSS_Solve(odeprob,mliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=100.0)
  #   solnliqss2=QSS_Solve(odeprob,nliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=100.0)
   solnmliqss2=QSS_Solve(odeprob,nmliqss1(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=100.0)



   #=  solliqss2Interp=solInterpolated(solliqss2,0.01)
    err1=getAverageErrorByRefs(solRodas5PVectorTyson,solliqss2Interp)
    #timeliqss=@belapsed QSS_Solve($odeprob,liqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resliqss= ("liqss",err1,solliqss2.totalSteps,solliqss2.simulStepCount,0.0)
    @show resliqss
    @btime QSS_Solve($odeprob,liqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)  =#
   #=  solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    err2=getAverageErrorByRefs(solRodas5PVectorTyson,solmliqss2Interp)
    timeliqss=@belapsed QSS_Solve($odeprob,liqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resmliqss= ("mliqss",err2,solmliqss2.totalSteps,solmliqss2.simulStepCount,timeliqss)
    @show resmliqss
=#
    #= solnliqss2Interp=solInterpolated(solnliqss2,0.01)
    err2=getAverageErrorByRefs(solRodas5PVectorTyson,solnliqss2Interp)
    timenliqss=@belapsed QSS_Solve($odeprob,nliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resnliqss= ("nliqss",err2,solnliqss2.totalSteps,solnliqss2.simulStepCount,timenliqss)
    @show resnliqss =#
    timenmliqss=0.0
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    err3=getAverageErrorByRefs(solRodas5PVectorTyson,solnmliqss2Interp)
    #timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
    resnmliqss= ("nmliqss",err3,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,timenmliqss)
    @show resnmliqss

  #=   XLSX.openxlsx("Tyson all solvers_$(relTol)_$(absTol).xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "Tyson all solvers__$(relTol)_$(absTol)"
      sheet["A4"] = collect(("solver","error","totalSteps","simul_steps","time"))
      sheet["A5"] = collect(resliqss)
     # sheet["A6"] = collect(resmliqss)
      sheet["A7"] = collect(resnliqss)
      sheet["A8"] = collect(resnmliqss)
   end =#
  
end

test53()
  