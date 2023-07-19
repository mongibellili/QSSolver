using formalqss
using BSON
using BenchmarkTools
using XLSX

function test()
 BSON.@load "formalqss/ref_bson/solVectAdvection_N10d1_Feagin14e-12.bson" solFeagin14VectorN10d1
  prob=@NLodeProblem begin
    name=(adrN10d01,)
    u[1:3]=1.0
    u[4:10]=0.0
    _dx=1.0#1/dx=N/10=10/10
    a=1.0
    d=0.1
    r=1000.0
   
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:9  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[10]=-a*_dx*(u[10]-u[9])+d*_dx*_dx*(2.0*u[9]-2.0*u[10])+r*u[10]*u[10]*(1.0-u[10]) 
    
end
    println("start solving")
   
 #=   solliqss2=QSS_Solve(prob,liqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    solmliqss2=QSS_Solve(prob,mliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
   solnliqss2=QSS_Solve(prob,nliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
   solnmliqss2=QSS_Solve(prob,nmliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    
    solliqss2Interp=solInterpolated(solliqss2,0.01)
    err1=getAverageErrorByRefs(solFeagin14VectorN10d1,solliqss2Interp)
    resliqss= ("liqss",err1,solliqss2.totalSteps,solliqss2.simulStepCount)
    @show resliqss
    
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    err=getAverageErrorByRefs(solFeagin14VectorN10d1,solmliqss2Interp)
    resmliqss= ("mliqss",err,solmliqss2.totalSteps,solmliqss2.simulStepCount)
    @show resmliqss

    solnliqss2Interp=solInterpolated(solnliqss2,0.01)
    err3=getAverageErrorByRefs(solFeagin14VectorN10d1,solnliqss2Interp)
    resnliqss= ("nliqss",err3,solnliqss2.totalSteps,solnliqss2.simulStepCount)
    @show resnliqss

    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    err=getAverageErrorByRefs(solFeagin14VectorN10d1,solnmliqss2Interp)
    resnmliqss= ("nmliqss",err,solnmliqss2.totalSteps,solnmliqss2.simulStepCount)
    @show resnmliqss


    XLSX.openxlsx("sysadrN10d1 all solvers______.xlsx", mode="w") do xf
        sheet = xf[1]
        sheet["A1"] = "adrN10d1all solvers_"
        sheet["A4"] = collect(("solver","error","totalSteps","simul_steps"))
        sheet["A5"] = collect(resliqss)
        sheet["A6"] = collect(resmliqss)
        sheet["A7"] = collect(resnliqss)
        sheet["A8"] = collect(resnmliqss)
     end


   



    
   # save_Sol(solnmliqss2)

   #@show length(solnmliqss2.savedTimes[1])

   
   for k=1:10
    save_Sol(solliqss2,"x$k",k;xlims=(0.0,10.0),ylims=(0.0,1.118))
    save_Sol(solmliqss2,"x$k",k;xlims=(0.0,10.0),ylims=(0.0,1.118))
    save_Sol(solnliqss2,"x$k",k;xlims=(0.0,10.0),ylims=(0.0,1.118))
    save_Sol(solnmliqss2,"x$k",k;xlims=(0.0,10.0),ylims=(0.0,1.118))
   # @show length(solnmliqss2.savedTimes[k])
   end =#


  # @btime QSS_Solve($prob,liqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)
  # @btime QSS_Solve($prob,mliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)
  # @btime QSS_Solve($prob,nliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)
  # @btime QSS_Solve($prob,nmliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)






#= 
   solqss2=QSS_Solve(prob,qss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    
   solliqss2Interp=solInterpolated(solqss2,0.01)
   err1=getAverageErrorByRefs(solFeagin14VectorN10d1,solliqss2Interp)
   resliqss= ("liqss",err1,solqss2.totalSteps,solqss2.simulStepCount)
   @show resliqss
   @btime QSS_Solve($prob,qss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0) =#

end
test()
