using formalqss
using BSON
using BenchmarkTools
using XLSX
function test()
    BSON.@load "formalqss/ref_bson/solVectAdvection_N100d1_Feagin14e-12.bson" solFeagin14VectorN100d1
  prob=@NLodeProblem begin
    name=(adrN100d1,)
    u[1:33]=1.0
    u[34:100]=0.0
    _dx=10.0#1/dx=N/10=100/10
    a=1.0
    d=1.0
    r=1000.0
   
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:99  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[100]=-a*_dx*(u[100]-u[99])+d*_dx*_dx*(2.0*u[99]-2.0*u[100])+r*u[100]*u[100]*(1.0-u[100]) 
    
end
    println("start solving")


    #solliqss2=QSS_Solve(prob,liqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
   #=  solliqss2Interp=solInterpolated(solliqss2,0.01)
    err1=getAverageErrorByRefs(solFeagin14VectorN100d1,solliqss2Interp)
    ttliqss=@belapsed QSS_Solve($prob,liqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)
    resliqss= ("liqss",err1,solliqss2.totalSteps,solliqss2.simulStepCount,ttliqss)
    @show resliqss =#

  #=   solmliqss2=QSS_Solve(prob,mliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    err2=getAverageErrorByRefs(solFeagin14VectorN100d1,solmliqss2Interp)
    ttmliqss=@belapsed QSS_Solve($prob,mliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)
    resmliqss= ("mliqss",err2,solmliqss2.totalSteps,solmliqss2.simulStepCount,ttmliqss)
    @show resmliqss

    solnliqss2=QSS_Solve(prob,nliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    solnliqss2Interp=solInterpolated(solnliqss2,0.01)
    err3=getAverageErrorByRefs(solFeagin14VectorN100d1,solnliqss2Interp)
    ttnliqss=@belapsed QSS_Solve($prob,nliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)
    resnliqss= ("nliqss",err3,solnliqss2.totalSteps,solnliqss2.simulStepCount,ttnliqss)
    @show resnliqss
=#
    solnmliqss2=QSS_Solve(prob,nmliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN100d1,solnmliqss2Interp)
    ttnmliqss=@belapsed QSS_Solve($prob,nmliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)
    resnmliqss= ("nmliqss",err4,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,ttnmliqss)
    @show resnmliqss

  
    #=
  
   
   
    XLSX.openxlsx("ADR N100d1.xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"] = "ADR_N100d1 all solvers_"
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
      
   end 
    =#
   
#=
    @btime solmliqss2=QSS_Solve($prob,liqss2()#= ,sparsity=sparse(),lightness=light() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    @btime solmliqss2=QSS_Solve($prob,mliqss2()#= ,sparsity=sparse(),lightness=light() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    @btime solmliqss2=QSS_Solve($prob,nliqss2()#= ,sparsity=sparse(),lightness=light() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    @btime solmliqss2=QSS_Solve($prob,nmliqss2()#= ,sparsity=sparse(),lightness=light() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
 =#




 
 #= solqss2=QSS_Solve(prob,qss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    
 solliqss2Interp=solInterpolated(solqss2,0.01)
 err1=getAverageErrorByRefs(solFeagin14VectorN100d1,solliqss2Interp)
 resliqss= ("liqss",err1,solqss2.totalSteps,solqss2.simulStepCount)
 @show resliqss
 @btime QSS_Solve($prob,qss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0) =#
end
test()
