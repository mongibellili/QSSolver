
using qss
using XLSX
#using BenchmarkTools
#using BSON

 function nl(absTol,relTol,solvr,opts,cycleDetection)
    prob=@NLodeProblem begin
      name=(nonlinear,)
      u  = [10.0,10.0,1.0,1.0]
      du[1] =  -0.1*u[1]*u[1]*u[2]+u[2]+1.0*u[4]
      du[2] =  -0.05*u[2]*u[2]*u[1]+u[1]-1.0*u[3]
      du[3] =  -0.05*u[3]+0.1
      du[4] =  -0.05*u[4]+0.1
    end
    #println("start solving")
    timenmliqss=0.0;err3=0.0
      timenmliqss=0.0
      tspan=(0.0,100.0)
       #  tspan=(0.0,0.17849696585)#58725 ---------
      sol=solve(prob,solvr,detection=Val(cycleDetection),options=opts,abstol=absTol,saveat=0.01,reltol=relTol,tspan,maxiters=20000000#= ,maxErr=100*relTol =#)
   #save_Sol(sol)
      #=   solInterp=solInterpolated(sol,0.01)
      err3=getAverageErrorByRefs(solInterp,solFeagin14VectorN1000d01) =#
     # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0,maxErr=10*$relTol)
      resnmliqss1E_2= ("$(sol.algName)",err3,sol.stats.totalSteps,sol.stats.simulStepCount,timenmliqss)
    # @show resnmliqss1E_2
      return resnmliqss1E_2
end 






function main(solver,multiplier,absTol,relTol)
         
  #BSON.@load "qss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
  
  solverName=(typeof(solver.name).parameters[1])
  note="NL_$(solverName)_abs$(absTol)_$(multiplier)Deltas"   # name of output file
  Results=[]
  for du=1:2
    for su=1:3
      for cycleDetection=0:4
          push!(Results,nl(absTol,relTol,solver,option(su,du,multiplier),cycleDetection))
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

multiplier=2.0
absTol=1e-4
relTol=1e-2
solver=mliqss1()
#nl(absTol,relTol,solver,option(1,1,multiplier),1)
main(solver,multiplier,absTol,relTol)
solver=nmliqss1()
#nl(absTol,relTol,solver,option(1,1,multiplier),1)
main(solver,multiplier,absTol,relTol)



#= multipliers=[1.0,2.0]
absTols=[1e-2,1e-3,1e-4]
relTols=1e-2

for multiplier in multipliers
  for absTol in absTols
    solver=mliqss1()
    main(solver,multiplier,absTol,relTol)
    solver=nmliqss1()
    main(solver,multiplier,absTol,relTol)
  end
end =#
