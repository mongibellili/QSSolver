
using qss
#using XLSX
#using BenchmarkTools
using BSON


function adr(case,slvr)
   # BSON.@load "OldQSS/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
    prob=@NLodeProblem begin
        name=(adrN1000d01,)
        u[1:333]=1.0
        u[334:1000]=0.0
        _dx=100.0#1/dx=N/10=1000/10
        a=1.0
        d=0.1
        r=1000.0
    
        du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
        for k in 2:999  
            du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
        end 
        du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
    end
    println("start solving")
    ttnmliqss=0.0

    absTol=1e-3
    relTol=1e-2
    sol=QSS_Solve(prob,slvr,dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
  #  solInterp=solInterpolated(sol,0.01)
 # err4=0.0
  #  err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solInterp)
    #ttnmliqss=@belapsed QSS_Solve($prob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
  #  resnmliqss1E_2= ("$(sol.algName)",relTol,err4,sol.totalSteps,sol.simulStepCount,ttnmliqss)
  #  @show resnmliqss1E_2
   #=  sol=QSS_Solve(prob,nmliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solInterp=solInterpolated(sol,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solInterp)
   # ttnmliqss=@belapsed QSS_Solve($prob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnmliqss2E_2= ("$(sol.algName)",relTol,err4,sol.totalSteps,sol.simulStepCount,ttnmliqss)
    @show resnmliqss2E_2 =#
   
 #=    absTol=1e-6
    relTol=1e-3
    sol=QSS_Solve(prob,nmliqss1(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solInterp=solInterpolated(sol,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solInterp)
    #ttnmliqss=@belapsed QSS_Solve($prob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnmliqss1E_3= ("$(sol.algName)",relTol,err4,sol.totalSteps,sol.simulStepCount,ttnmliqss)
    @show resnmliqss1E_3
    sol=QSS_Solve(prob,nmliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solInterp=solInterpolated(sol,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solInterp)
   # ttnmliqss=@belapsed QSS_Solve($prob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnmliqss2E_3= ("$(sol.algName)",relTol,err4,sol.totalSteps,sol.simulStepCount,ttnmliqss)
    @show resnmliqss2E_3

    absTol=1e-7
    relTol=1e-4
    sol=QSS_Solve(prob,nmliqss1(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solInterp=solInterpolated(sol,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solInterp)
    #ttnmliqss=@belapsed QSS_Solve($prob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnmliqss1E_4= ("$(sol.algName)",relTol,err4,sol.totalSteps,sol.simulStepCount,ttnmliqss)
    @show resnmliqss1E_4
    sol=QSS_Solve(prob,nmliqss2(),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
    solInterp=solInterpolated(sol,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solInterp)
   # ttnmliqss=@belapsed QSS_Solve($prob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnmliqss2E_4= ("$(sol.algName)",relTol,err4,sol.totalSteps,sol.simulStepCount,ttnmliqss)
    @show resnmliqss2E_4
  
    XLSX.openxlsx("ADR N1000d01_$case.xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "ADR_N1000d01 $case"
      sheet["A4"] = collect(("solver","Tolerance","error","totalSteps","simul_steps","time"))
      sheet["A5"] = collect(resnmliqss1E_2)
      sheet["A6"] = collect(resnmliqss1E_3)
      sheet["A7"] = collect(resnmliqss1E_4)

      sheet["A8"] = collect(resnmliqss2E_2)
      sheet["A9"] = collect(resnmliqss2E_3)
      sheet["A10"] = collect(resnmliqss2E_4)
    end
  =#
  
    
   
end 

case="bothMaxTol"
#test53(case)
#tyson(case)
@time adr(case,mliqss1())
@time adr(case,nmliqss1())