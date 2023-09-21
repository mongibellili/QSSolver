
using qss
#using XLSX
#using BenchmarkTools
using BSON
#using TimerOutputs
#using Plots
function test(case,solvr)
  absTol=1e-5
     relTol=1e-2
      odeprob = @NLodeProblem begin
         #sys b53
         name=(sysb53,)
         u = [-1.0, -2.0]
         du[1] = -20.0*u[1]-80.0*u[2]+1600.0
         du[2] =1.24*u[1]-0.01*u[2]+0.2
     end  
     u1, u2 = -8.73522174738572, -7.385745994549763
     λ1, λ2 = -10.841674966758294, -9.168325033241706
     c1, c2 = 121.14809142478035, -143.14809142478035
     xp1, xp2 = 0.0, 20.0
     x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
     x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

     timenmliqss=0.0;er1=0.0;er2=0.0
     
     println("start LTI solving")
     solnmliqss=QSS_Solve(odeprob,solvr,dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=100.0,maxErr=10000*relTol)

     solnmliqssInterp=solInterpolated(solnmliqss,0.01)
     er1=getError(solnmliqssInterp,1,x1)  
     er2=getError(solnmliqssInterp,2,x2) 
   # timenmliqss=@belapsed QSS_Solve($odeprob,$solvr,dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0#= ,maxErr=1000*$relTol =#)
     resnmliqss1E_2= ("$(solnmliqss.algName)",relTol,(er1+er2)/2,solnmliqss.totalSteps,solnmliqss.simulStepCount,timenmliqss)
     @show resnmliqss1E_2



    
      BSON.@load "qss/ref_bson/solVect_Tyson_Rodas5Pe-12.bson" solRodas5PVectorTyson
     odeprob = @NLodeProblem begin
         name=(tyson,)
         u = [0.0,0.75,0.25,0.0,0.0,0.0]
         du[1] = u[4]-1e6*u[1]+1e3*u[2]
         du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
         du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
         du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
         du[5] = 0.015-200.0*u[2]*u[5]
         du[6] =u[4]-0.6*u[6]
     end   
     println("start tyson solving")
     timenmliqss=0.0
     solnmliqss=QSS_Solve(odeprob,solvr,dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0#= ,maxErr=100*relTol =#)
    # save_Sol(solnmliqss,1,note="x1 intrval13  ",xlims=(4.0,4.38),ylims=(0.0007,0.000723))
    solnmliqssInterp=solInterpolated(solnmliqss,0.01)
    err3=getAverageErrorByRefs(solRodas5PVectorTyson,solnmliqssInterp) 

  # timenmliqss=@belapsed QSS_Solve($odeprob,$solvr,dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0#= ,maxErr=1000*$relTol =#)
    resnmliqss11E_2= ("$(solnmliqss.algName)",relTol,err3,solnmliqss.totalSteps,solnmliqss.simulStepCount,timenmliqss)
    @show resnmliqss11E_2 
     


      BSON.@load "qss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
    prob=@NLodeProblem begin
        name=(adrN1000d01,)
        u[1:333]=1.0
        u[334:1000]=0.0
        _dx=100.0#1/dx=N/10=1000/10
        a=1.0;d=0.1;r=1000.0
        du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
        for k in 2:999  
            du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
        end 
        du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
    end
    println("start adr solving")
    ttnmliqss=0.0

  
    solnmliqss=QSS_Solve(prob,solvr,dQmin=absTol,saveat=0.01,dQrel=relTol,finalTime=10.0)#
   # save_Sol(solnmliqss,1,2,200,400,600)
     solnmliqssInterp=solInterpolated(solnmliqss,0.01)
    err4=getAverageErrorByRefs(solFeagin14VectorN1000d01,solnmliqssInterp)
   #ttnmliqss=@belapsed QSS_Solve($prob,$solvr,dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=10.0)
    resnmliqss12E_2= ("$(solnmliqss.algName)",relTol,err4,solnmliqss.totalSteps,solnmliqss.simulStepCount,ttnmliqss)
    @show resnmliqss12E_2    



    #=  XLSX.openxlsx("3sys $(solvr)_$(case)_$(relTol).xlsx", mode="w") do xf
        sheet = xf[1]
        sheet["A1"] = "3sys __$case)"
        sheet["A4"] = collect(("solver","Tolerance","error","totalSteps","simul_steps","time"))
        sheet["A5"] = collect(resnmliqss1E_2)
        sheet["A6"] = collect(resnmliqss11E_2)
        sheet["A7"] = collect(resnmliqss12E_2)

       #=  sheet["A8"] = collect(resnmliqss2E_2)
        sheet["A9"] = collect(resnmliqss2E_3)
        sheet["A10"] = collect(resnmliqss2E_4) =#
     end =#
end

case="order1_newcondsijdxithrow"
#test(case,mliqss1())
test(case,nmliqss1())