
using qss
using XLSX
using BenchmarkTools
using BSON
function adr(absTol,relTol,solFeagin14VectorN1000d01,solvr,opts,cycleDetection)
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
    timenmliqss=0.0;err3=0.0
    tspan=(0.0,10.0)
    sol=solve(prob,solvr,detection=Val(cycleDetection),options=opts,abstol=absTol,saveat=0.01,reltol=relTol,tspan,maxiters=20000000#= ,maxErr=100*relTol =#)
    #save_Sol(sol,1)
    solInterp=solInterpolated(sol,0.01)
    err3=getAverageErrorByRefs(solInterp,solFeagin14VectorN1000d01)
    timenmliqss=@belapsed  solve($prob,$solvr,detection=Val($cycleDetection),options=$opts,abstol=$absTol,reltol=$relTol,$tspan,maxiters=20000000#= ,maxErr=100*relTol =#)
 
    resnmliqss1E_2= ("$(sol.algName)",err3,sol.stats.totalSteps,sol.stats.simulStepCount,timenmliqss)
    @show resnmliqss1E_2
    return resnmliqss1E_2  
end 

function main(solFeagin14VectorN1000d01,solver,multiplier,absTol,relTol)
  println("__________solver=$solver")
  #BSON.@load "qss/ref_bson/solVect_adr_Rodas5Pe-12" solFeagin14VectorN1000d01
  
  solverName=(typeof(solver.name).parameters[1])
  note="adr_$(solverName)_abs$(absTol)_$(multiplier)_run2_singles_OldCoeff"   # name of output file
  Results=[]
  for du=1:1 # du double update 1 iters, 2 analytical
    for su=3:3 #su single update 1,2 analytical, 3 iters
      for cycleDetection=1:4 #0,1,2,3,4
          push!(Results,adr(absTol,relTol,solFeagin14VectorN1000d01,solver,option(su,du,multiplier),cycleDetection))
      end
    end
  end
 #=  push!(Results,adr(absTol,relTol,solFeagin14VectorN1000d01,solver,option(1,1,multiplier),0))
  push!(Results,adr(absTol,relTol,solFeagin14VectorN1000d01,solver,option(1,2,multiplier),2))
  push!(Results,adr(absTol,relTol,solFeagin14VectorN1000d01,solver,option(2,2,multiplier),4))
  push!(Results,adr(absTol,relTol,solFeagin14VectorN1000d01,solver,option(3,2,multiplier),2)) =#
#=   XLSX.openxlsx("_$(note)_.xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "_$(note)_"
    sheet["A2"] = collect(("solver and options","error","totalSteps","simul_steps","time"))
    for i=3:length(Results)+2
      sheet["A$i"] = collect(Results[i-2])
    end
  end =#


end 


BSON.@load "qss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01

multipliers=[1#= ,4,5,6,7,8,9 =#] # to test maxiters in while loop...not needed any more just use 1. 
absTols=[1e-2]
relTol=1e-1
for multiplier in multipliers
  for absTol in absTols
    println("____absTol=$absTol")
    solver=mliqss1() # mliqss is related to original calculation of linear coefficients
    main(solFeagin14VectorN1000d01,solver,multiplier,absTol,relTol)
    solver=nmliqss1() # nmliqss is related to new calculation of linear coefficients
    main(solFeagin14VectorN1000d01,solver,multiplier,absTol,relTol)
  end
end
