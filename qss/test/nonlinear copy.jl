
using qss
using XLSX
#using BenchmarkTools
using BSON

 function Oregon(solRodas5PVectorOregon,note,slvr,cycleDetection)
  
    
       odeprob = @NLodeProblem begin
            name=(nonlinear,)
            u  = [10.0,10.0]
          
            #= du[1] =  -0.05*u[1]
            du[2] =  -0.05*u[2] =#
            du[1] =  -0.1*u[1]*u[1]*u[2]+u[2]
            du[2] =  -0.05*u[2]*u[2]*u[1]+u[1]

            
          
       end  
     
    
       println("start solving")
       timenmliqss=0.0
     
       absTol=1e-3
       relTol=1e-2
     
      sol=QSS_Solve(odeprob,slvr,detection=Val(cycleDetection),dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=100.0#= ,maxErr=10*relTol =#)
    # save_Sol(sol)
      #solInterp=solInterpolated(sol,0.01)
      err3=0.0
      #err3=getAverageErrorByRefs(solRodas5PVectorOregon,solInterp)
     # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0,maxErr=10*$relTol)
      resnmliqss1E_2= ("$(sol.algName)",err3,sol.totalSteps,sol.simulStepCount,timenmliqss)
      @show resnmliqss1E_2

      return resnmliqss1E_2
    
end 
solRodas5PVectorOregon=0.0
note="abs-3_1Delta" 
resmliqss21=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulIter(),1)
resmliqss22=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulIter(),2)
#= resmliqss23=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulIter(),3)
resmliqss24=Oregon(solRodas5PVectorOregon,note,mliqss2_SimulIter(),4) =#
