using formalqss
using BSON
using BenchmarkTools
using XLSX

function test()
  BSON.@load "formalqss/ref_bson/solVectAdvection_N1000d01r10_Feagin14e-12.bson" solFeagin14VectorN1000d01r10  
  prob=@NLodeProblem begin
    u[1:333]=1.0
    u[334:1000]=0.0
    _dx=100.0#1/dx=N/10=1000/10
    a=1.0
    d=0.1
    r=10.0
   
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:999  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
    
end
    println("start solving")
    absTol=1e-5
    relTol=1e-2
    solliqss2=QSS_Solve(prob,liqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solliqss2Interp=solInterpolated(solliqss2,0.01)
    err1=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solliqss2Interp)
    ttliqss=@belapsed QSS_Solve($prob,liqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resliqss_1= ("liqss",err1,solliqss2.totalSteps,solliqss2.simulStepCount,ttliqss)
    @show resliqss_1

    solnliqss2=QSS_Solve(prob,nliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    @show solnliqss2.totalSteps
    solnliqss2Interp=solInterpolated(solnliqss2,0.01)
    err3=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solnliqss2Interp)
    ttnliqss=@belapsed QSS_Solve($prob,nliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnliqss_1= ("nliqss",err3,solnliqss2.totalSteps,solnliqss2.simulStepCount,ttnliqss)
    @show resnliqss_1
    
    solnmliqss2=QSS_Solve(prob,nmliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solnmliqss2Interp)
    ttnmliqss=@belapsed QSS_Solve($prob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnmliqss_1= ("nmliqss",err4,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,ttnmliqss)
    @show resnmliqss_1

  #=   solqss2=QSS_Solve(prob,qss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solqss2Interp=solInterpolated(solqss2,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solqss2Interp)
    ttqss=@belapsed QSS_Solve($prob,qss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resqss_1= ("qss",err4,solqss2.totalSteps,solqss2.simulStepCount,ttqss)
    @show resqss_1 =#
 
    absTol=1e-6
    relTol=1e-3
 solliqss2=QSS_Solve(prob,liqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
 solliqss2Interp=solInterpolated(solliqss2,0.01)
 err1=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solliqss2Interp)
 ttliqss=@belapsed QSS_Solve($prob,liqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
 resliqss_2= ("liqss",err1,solliqss2.totalSteps,solliqss2.simulStepCount,ttliqss)
 @show resliqss_2

 solnliqss2=QSS_Solve(prob,nliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
 @show solnliqss2.totalSteps
 solnliqss2Interp=solInterpolated(solnliqss2,0.01)
 err3=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solnliqss2Interp)
 ttnliqss=@belapsed QSS_Solve($prob,nliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
 resnliqss_2= ("nliqss",err3,solnliqss2.totalSteps,solnliqss2.simulStepCount,ttnliqss)
 @show resnliqss_2
 
 solnmliqss2=QSS_Solve(prob,nmliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
 solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
 err4=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solnmliqss2Interp)
 ttnmliqss=@belapsed QSS_Solve($prob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
 resnmliqss_2= ("nmliqss",err4,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,ttnmliqss)
 @show resnmliqss_2

 #= solqss2=QSS_Solve(prob,qss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
 solqss2Interp=solInterpolated(solqss2,0.01)
 err4=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solqss2Interp)
 ttqss=@belapsed QSS_Solve($prob,qss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
 resqss_2= ("qss",err4,solqss2.totalSteps,solqss2.simulStepCount,ttqss)
 @show resqss_2 =#
    

 absTol=1e-7
    relTol=1e-4
    solliqss2=QSS_Solve(prob,liqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solliqss2Interp=solInterpolated(solliqss2,0.01)
    err1=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solliqss2Interp)
    ttliqss=@belapsed QSS_Solve($prob,liqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resliqss_3= ("liqss",err1,solliqss2.totalSteps,solliqss2.simulStepCount,ttliqss)
    @show resliqss_3
   
    solnliqss2=QSS_Solve(prob,nliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    @show solnliqss2.totalSteps
    solnliqss2Interp=solInterpolated(solnliqss2,0.01)
    err3=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solnliqss2Interp)
    ttnliqss=@belapsed QSS_Solve($prob,nliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnliqss_3= ("nliqss",err3,solnliqss2.totalSteps,solnliqss2.simulStepCount,ttnliqss)
    @show resnliqss_3
    
    solnmliqss2=QSS_Solve(prob,nmliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solnmliqss2Interp)
    ttnmliqss=@belapsed QSS_Solve($prob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnmliqss_3= ("nmliqss",err4,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,ttnmliqss)
    @show resnmliqss_3
   
    #= solqss2=QSS_Solve(prob,qss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solqss2Interp=solInterpolated(solqss2,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01r10,solqss2Interp)
    ttqss=@belapsed QSS_Solve($prob,qss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resqss_3= ("qss",err4,solqss2.totalSteps,solqss2.simulStepCount,ttqss)
    @show resqss_3 =#
  
    XLSX.openxlsx("ADR N1000d01r10.xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "ADR_N1000d01r10 all solvers_"
      sheet["A2"] = "e-2_e-5_"
      sheet["A3"] = collect(("solver","error","totalSteps","simul_steps","time"))
      sheet["A4"] = collect(resliqss_1)
      sheet["A5"] = collect(resnliqss_1)
      sheet["A6"] = collect(resnmliqss_1)
      #sheet["A7"] = collect(resqss_1)
  
      sheet["A8"] = "e-3_e-6_"
      sheet["A9"] = collect(("solver","error","totalSteps","simul_steps","time"))
      sheet["A10"] = collect(resliqss_2)
      sheet["A11"] = collect(resnliqss_2)
      sheet["A12"] = collect(resnmliqss_2)
      #sheet["A13"] = collect(resqss_2)

     

      sheet["A14"] = "e-4_e-7_"
      sheet["A15"] = collect(("solver","error","totalSteps","simul_steps","time"))
      sheet["A16"] = collect(resliqss_3)
      sheet["A17"] = collect(resnliqss_3)
      sheet["A18"] = collect(resnmliqss_3)
      #sheet["A19"] = collect(resqss_3)

      #= for k=1:1000
        sheet["A$(3*k+4)"]="variable $k changes: $(length(solnmliqss2.savedTimes[k]))"
        sheet["A$(3*k+5)"]=collect(solnmliqss2.savedTimes[k])
        sheet["A$(3*k+6)"]=collect(solnmliqss2.savedVars[k])
      end =#
    end
end
test()
