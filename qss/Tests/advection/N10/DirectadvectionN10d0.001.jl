using formalqss
using BSON
using BenchmarkTools

function test()
  BSON.@load "/home/mongi/projects/formalqss/solVectAdvection_N10d001_Feagin14e-12.bson" solFeagin14VectorN10d001
  prob=@NLodeProblem begin
    name=(adrN10d001,)
    u[1:3]=1.0
    u[4:10]=0.0
    _dx=1.0#1/dx=N/10=10/10
    a=1.0
    d=0.001
    r=1000.0
   
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:9  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[10]=-a*_dx*(u[10]-u[9])+d*_dx*_dx*(2.0*u[9]-2.0*u[10])+r*u[10]*u[10]*(1.0-u[10]) 
    
end
    println("start solving")
   
    solliqss2=QSS_Solve(prob,liqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    solliqss2Interp=solInterpolated(solliqss2,0.01)
    err1=getAverageErrorByRefs(solFeagin14VectorN10d001,solliqss2Interp)
    resliqss= ("liqss",err1,solliqss2.totalSteps,solliqss2.simulStepCount)
    @show resliqss

    solmliqss2=QSS_Solve(prob,mliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    err=getAverageErrorByRefs(solFeagin14VectorN10d001,solmliqss2Interp)
    resmliqss= ("mliqss",err,solmliqss2.totalSteps,solmliqss2.simulStepCount)
    @show resmliqss

    solnliqss2=QSS_Solve(prob,nliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    solnliqss2Interp=solInterpolated(solnliqss2,0.01)
    err3=getAverageErrorByRefs(solFeagin14VectorN10d001,solnliqss2Interp)
    resnliqss= ("nliqss",err3,solnliqss2.totalSteps,solnliqss2.simulStepCount)
    @show resnliqss

    solmliqss2=QSS_Solve(prob,nmliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    err=getAverageErrorByRefs(solFeagin14VectorN10d001,solmliqss2Interp)
    resnmliqss= ("nmliqss",err,solmliqss2.totalSteps,solmliqss2.simulStepCount)
    @show resnmliqss

    @btime solmliqss2=QSS_Solve($prob,liqss2()#= ,sparsity=sparse(),lightness=light() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    @btime solmliqss2=QSS_Solve($prob,mliqss2()#= ,sparsity=sparse(),lightness=light() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    @btime solmliqss2=QSS_Solve($prob,nliqss2()#= ,sparsity=sparse(),lightness=light() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    @btime solmliqss2=QSS_Solve($prob,nmliqss2()#= ,sparsity=sparse(),lightness=light() =#,dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#

end
test()
