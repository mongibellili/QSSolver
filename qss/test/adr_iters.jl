
using qss
using XLSX
#using BenchmarkTools
using BSON

 function adr(solFeagin14VectorN1000d01,note,slvr,cycleDetection)
  
    
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
sol=QSS_Solve(prob,slvr,detection=Val(cycleDetection),dQmin=absTol,dQrel=relTol,finalTime=10.0)#
  solInterp=solInterpolated(sol,0.01)
 err4=0.0
  err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solInterp)
#ttnmliqss=@belapsed QSS_Solve($prob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
  resnmliqss1E_2= ("$(sol.algName)",err4,sol.totalSteps,sol.simulStepCount,ttnmliqss)
  @show resnmliqss1E_2

      return resnmliqss1E_2
    
end 

BSON.@load "qss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
note="iters_abs-3"
#= resmliqss10=adr(solFeagin14VectorN1000d01,note,mliqss1_SimulIter(),0)
resmliqss11=adr(solFeagin14VectorN1000d01,note,mliqss1_SimulIter(),1)
resmliqss12=adr(solFeagin14VectorN1000d01,note,mliqss1_SimulIter(),2)
resmliqss13=adr(solFeagin14VectorN1000d01,note,mliqss1_SimulIter(),3)
resmliqss14=adr(solFeagin14VectorN1000d01,note,mliqss1_SimulIter(),4) =#

resmliqss20=adr(solFeagin14VectorN1000d01,note,mliqss2_SimulIter(),0)
resmliqss21=adr(solFeagin14VectorN1000d01,note,mliqss2_SimulIter(),1)
resmliqss22=adr(solFeagin14VectorN1000d01,note,mliqss2_SimulIter(),2)
resmliqss23=adr(solFeagin14VectorN1000d01,note,mliqss2_SimulIter(),3)
resmliqss24=adr(solFeagin14VectorN1000d01,note,mliqss2_SimulIter(),4)

resmliqss30=adr(solFeagin14VectorN1000d01,note,mliqss3_SimulIter(),0)
resmliqss31=adr(solFeagin14VectorN1000d01,note,mliqss3_SimulIter(),1)
resmliqss32=adr(solFeagin14VectorN1000d01,note,mliqss3_SimulIter(),2)
resmliqss33=adr(solFeagin14VectorN1000d01,note,mliqss3_SimulIter(),3)
resmliqss34=adr(solFeagin14VectorN1000d01,note,mliqss3_SimulIter(),4)

XLSX.openxlsx("adr_solvers_$note.xlsx", mode="w") do xf
  sheet = xf[1]
  sheet["A1"] = "adr_solvers__$note)"
  sheet["A4"] = collect(("solver","error","totalSteps","simul_steps","time"))
#=   sheet["A5"] = collect(resmliqss10)
  sheet["A6"] = collect(resmliqss11)
  sheet["A7"] = collect(resmliqss12)
  sheet["A8"] = collect(resmliqss13)
  sheet["A9"] = collect(resmliqss14) =#

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