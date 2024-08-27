
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
    #println("start solving")
    timenmliqss=0.0;err3=0.0
     
      
     
     
      timenmliqss=0.0
      tspan=(0.0,25.0)
       #  tspan=(0.0,0.17849696585)#58725 ---------
     
  
      sol=solve(prob,solvr,detection=Val(cycleDetection),options=opts,abstol=absTol,saveat=0.01,reltol=relTol,tspan,maxiters=20000000#= ,maxErr=100*relTol =#)
  #save_Sol(sol,1,ylims=(0.0,0.01))
   #    solInterp=solInterpolated(sol,0.01)
    #  err3=getAverageErrorByRefs(solInterp,solRodas5PVectorTyson)
     # timenmliqss=@belapsed  solve($prob,$solvr,detection=Val($cycleDetection),options=$opts,abstol=$absTol,reltol=$relTol,$tspan,maxiters=20000000#= ,maxErr=100*relTol =#)
 
      resnmliqss1E_2= ("$(sol.algName)",err3,sol.stats.totalSteps,sol.stats.simulStepCount,timenmliqss)
     # @show resnmliqss1E_2

      return resnmliqss1E_2
    
    
end 






function main(solRodas5PVectorTyson,solver,multiplier,absTol,relTol)
  println("__________solver=$solver")
  #BSON.@load "qss/ref_bson/solVect_Tyson_Rodas5Pe-12" solRodas5PVectorTyson
  
  solverName=(typeof(solver.name).parameters[1])
  note="tyson_$(solverName)_abs$(absTol)_$(multiplier)Deltas_du2_detect1213"   # name of output file
  Results=[]
  for du=1:2
    for su=1:3
      for cycleDetection=0:4
          push!(Results,tyson(absTol,relTol,solRodas5PVectorTyson,solver,option(su,du,multiplier),cycleDetection))
      end
    end
  end
 #=  push!(Results,tyson(absTol,relTol,solRodas5PVectorTyson,solver,option(1,1,multiplier),0))
  push!(Results,tyson(absTol,relTol,solRodas5PVectorTyson,solver,option(1,2,multiplier),2))
  push!(Results,tyson(absTol,relTol,solRodas5PVectorTyson,solver,option(2,2,multiplier),4))
  push!(Results,tyson(absTol,relTol,solRodas5PVectorTyson,solver,option(3,2,multiplier),2)) =#
  XLSX.openxlsx("_$(note)_.xlsx", mode="w") do xf
    sheet = xf[1]
    sheet["A1"] = "_$(note)_"
    sheet["A2"] = collect(("solver and options","error","totalSteps","simul_steps","time"))
    for i=3:length(Results)+2
      sheet["A$i"] = collect(Results[i-2])
    end
  end


end 


BSON.@load "qss/ref_bson/solVect_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson


#= 
multiplier=2.0
absTol=1e-4
relTol=1e-2
solver=mliqss1()
#tyson(absTol,relTol,solRodas5PVectorTyson,solver,option(1,1,multiplier),1)
#main(solRodas5PVectorTyson,solver,multiplier,absTol,relTol)
solver=nmliqss1()
tyson(absTol,relTol,solRodas5PVectorTyson,solver,option(3,2,multiplier),4)
#main(solRodas5PVectorTyson,solver,multiplier,absTol,relTol) =#



multipliers=[2.0]
#= absTols=[1e-3,1e-4]
relTol=1e-2 =#
absTols=[1e-3]
relTol=1e-2

for multiplier in multipliers
  for absTol in absTols
    println("____absTol=$absTol")
    solver=mliqss1()
    main(solRodas5PVectorTyson,solver,multiplier,absTol,relTol)
    solver=nmliqss1()
    main(solRodas5PVectorTyson,solver,multiplier,absTol,relTol)
  end
end
