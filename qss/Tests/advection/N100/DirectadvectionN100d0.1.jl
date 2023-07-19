using formalqss
using BSON
using BenchmarkTools
using XLSX
function test()
    BSON.@load "formalqss/ref_bson/solVectAdvection_N100d01_Feagin14e-12.bson" solFeagin14VectorN100d01
  prob=@NLodeProblem begin
    name=(adrN100d01,)
    u[1:33]=1.0
    u[34:100]=0.0
    _dx=10.0#1/dx=N/10=100/10
    a=1.0
    d=0.1
    r=1000.0
   
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:99  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[100]=-a*_dx*(u[100]-u[99])+d*_dx*_dx*(2.0*u[99]-2.0*u[100])+r*u[100]*u[100]*(1.0-u[100]) 
    
end
    println("start solving")

ft=10.0
   #=  solliqss2=QSS_Solve(prob,liqss3(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft)#
    @show solliqss2.totalSteps 
    solliqss2Interp=solInterpolated(solliqss2,0.01)
    err1=getAverageErrorByRefs(solFeagin14VectorN100d01,solliqss2Interp)
   # ttliqss=@belapsed QSS_Solve($prob,liqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft)
    resliqss= ("liqss",err1,solliqss2.totalSteps,solliqss2.simulStepCount,0.0)
    @show resliqss =#

   #=  solmliqss2=QSS_Solve(prob,mliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft)#
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    err2=getAverageErrorByRefs(solFeagin14VectorN100d01,solmliqss2Interp)
    ttmliqss=@belapsed QSS_Solve($prob,mliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft)
    resmliqss= ("mliqss",err2,solmliqss2.totalSteps,solmliqss2.simulStepCount,ttmliqss)
    @show resmliqss

    solnliqss2=QSS_Solve(prob,nliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft)#
    solnliqss2Interp=solInterpolated(solnliqss2,0.01)
    err3=getAverageErrorByRefs(solFeagin14VectorN100d01,solnliqss2Interp)
    ttnliqss=@belapsed QSS_Solve($prob,nliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft)
    resnliqss= ("nliqss",err3,solnliqss2.totalSteps,solnliqss2.simulStepCount,ttnliqss)
    @show resnliqss =#

    solnmliqss2=QSS_Solve(prob,nmliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft)#
   @show  solnmliqss2.totalSteps,solnmliqss2.simulStepCount


  # save_Sol(solnmliqss2,"x333",33,34,35;xlims=(0.0,100.0),ylims=(0.0,1.018))


  # save_Sol(solliqss2,"x333",33,34,35;xlims=(0.0,10.0),ylims=(0.0,1.018))
     solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN100d01,solnmliqss2Interp)
   # ttnmliqss=@belapsed QSS_Solve($prob,nmliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft)
    resnmliqss= ("nmliqss",err4,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,0.0)
    @show resnmliqss

  
  
   
   
   #=  XLSX.openxlsx("ADR N100d01.xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "ADR_N100d01 all solvers_"
      sheet["A2"] = collect(("solver","error","totalSteps","simul_steps","time"))
      sheet["A3"] = collect(resliqss)
      sheet["A4"] = collect(resmliqss)
      sheet["A5"] = collect(resnliqss)
      sheet["A6"] = collect(resnmliqss)
  
      for k=1:100
        sheet["A$(3*k+4)"]="variable $k changes: $(length(solnmliqss2.savedTimes[k]))"
        sheet["A$(3*k+5)"]=collect(solnmliqss2.savedTimes[k])
        sheet["A$(3*k+6)"]=collect(solnmliqss2.savedVars[k])
      end
      
   end  =#
   
   

   #=  @btime solmliqss2=QSS_Solve($prob,liqss3()#= ,sparsity=sparse() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=$ft)#
 #   @btime solmliqss2=QSS_Solve($prob,mliqss2(),sparsity=sparse(),lightness=light(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=$ft)#
   # @btime solmliqss2=QSS_Solve($prob,nliqss3()#= ,sparsity=sparse() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=$ft)#
    @btime solmliqss2=QSS_Solve($prob,nmliqss3()#= ,sparsity=sparse() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=$ft)# =#





 
 #= solqss2=QSS_Solve(prob,qss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft)#
    
 solliqss2Interp=solInterpolated(solqss2,0.01)
 err1=getAverageErrorByRefs(solFeagin14VectorN100d01,solliqss2Interp)
 resliqss= ("liqss",err1,solqss2.totalSteps,solqss2.simulStepCount)
 @show resliqss
 @btime QSS_Solve($prob,qss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=ft) =#
end
test()
