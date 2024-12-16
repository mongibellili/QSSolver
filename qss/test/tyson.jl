
using qss
using XLSX
using BenchmarkTools
using BSON
function tyson(absTol,relTol,solRodas5PVectorTyson,solvr,opts,cycleDetection)
  prob=@NLodeProblem begin
    name=(tyson,)
           u = [0.0,0.75,0.25,0.0,0.0,0.0]
           du[1] = u[4]-1e6*u[1]+1e3*u[2]
           du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
           du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
           du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
           du[5] = 0.015-200.0*u[2]*u[5]
           du[6] =u[4]-0.6*u[6]
    end
    timenmliqss=0.0;err3=0.0
    tspan=(0.0,25.0)  
    sol=solve(prob,solvr,detection=Val(cycleDetection),options=opts,abstol=absTol,reltol=relTol,tspan,maxiters=20000000#= ,maxErr=100*relTol =#)
     # save_Sol(sol,1)
     # save_Sol(sol,5)
      #save_Sol(sol,1,ylims=(0.0,0.01))
    solInterp=solInterpolated(sol,0.01)
    err3=getAverageErrorByRefs(solInterp,solRodas5PVectorTyson)
     # timenmliqss=@belapsed  solve($prob,$solvr,detection=Val($cycleDetection),options=$opts,abstol=$absTol,reltol=$relTol,$tspan,maxiters=20000000#= ,maxErr=100*relTol =#)
    resnmliqss1E_2= ("$(sol.algName)",err3,sol.stats.totalSteps,sol.stats.simulStepCount,timenmliqss)
    @show resnmliqss1E_2
    return resnmliqss1E_2
end 

function main(solRodas5PVectorTyson,solver,multiplier,absTol,relTol)
  println("__________solver=$solver")
  solverName=(typeof(solver.name).parameters[1])
  note="tyson_$(solverName)_abs$(absTol)_$(multiplier)_summary"   # name of output file
  Results=[]
  for du=1:1
    for su=1:1
      for cycleDetection=1:1
          push!(Results,tyson(absTol,relTol,solRodas5PVectorTyson,solver,option(su,du,multiplier),cycleDetection))
      end
    end
  end

 #=  XLSX.openxlsx("_$(note)_.xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "_$(note)_"
    sheet["A2"] = collect(("solver and options","error","totalSteps","simul_steps","time"))
    for i=3:length(Results)+2
      sheet["A$i"] = collect(Results[i-2])
    end
  end =#
end 


BSON.@load "qss/ref_bson/solRodas5PVectorTyson.bson" solRodas5PVectorTyson
multipliers=[1]
absTols=[1e-5]
relTol=1e-2
for multiplier in multipliers
  for absTol in absTols
    println("____absTol=$absTol")
    solver=mliqss1()
    main(solRodas5PVectorTyson,solver,multiplier,absTol,relTol)
  #=   solver=nmliqss1()
    main(solRodas5PVectorTyson,solver,multiplier,absTol,relTol) =#
  end
end
