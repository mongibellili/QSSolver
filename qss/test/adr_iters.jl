
using qss
using XLSX
#using BenchmarkTools
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
println("start solving")
timenmliqss=0.0;err3=0.0
     
      
     
     
      timenmliqss=0.0
      tspan=(0.0,10.0)
       #  tspan=(0.0,0.17849696585)#58725 ---------
     
  
      sol=solve(prob,solvr,detection=Val(cycleDetection),options=opts,abstol=absTol,saveat=0.01,reltol=relTol,tspan,maxiters=20000000#= ,maxErr=100*relTol =#)
      solInterp=solInterpolated(sol,0.01)
      err3=getAverageErrorByRefs(solInterp,solFeagin14VectorN1000d01)
     # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0,maxErr=10*$relTol)
      resnmliqss1E_2= ("$(sol.algName)",err3,sol.stats.totalSteps,sol.stats.simulStepCount,timenmliqss)
      @show resnmliqss1E_2

      return resnmliqss1E_2
    
    
end 






function main(solver,multiplier,absTol,relTol)
         
  BSON.@load "qss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
  
  solverName=(typeof(solver.name).parameters[1])
  note="_$(solverName)_abs$(absTol)_$(multiplier)Deltas"   # name of output file
  Results=[]
  for du=1:1
    for su=1:1
      for cycleDetection=0:1
          push!(Results,adr(absTol,relTol,solFeagin14VectorN1000d01,solver,option(su,du,multiplier),cycleDetection))
      end
    end
  end
  XLSX.openxlsx("_$(note)_.xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "_$(note)_"
    sheet["A2"] = collect(("solver and options","error","totalSteps","simul_steps","time"))
    for i=3:length(Results)+2
      sheet["A$i"] = collect(Results[i-2])
    end
  end
end

multipliers=[1.0,2.0]
absTols=[1e-2,1e-3,1e-4]
relTols=1e-2

for multiplier in multipliers
  for absTol in absTols
    solver=mliqss1()
    main(solver,multiplier,absTol,relTol)
    solver=nmliqss1()
    main(solver,multiplier,absTol,relTol)
  end
end
