
using qss
using XLSX
#using BenchmarkTools

function test53(case,slvr,cycleDetection)
     odeprob = @NLodeProblem begin
         #sys b53
         name=(sysb53,)
         u = [-1.0, -2.0]
        # discrete = [0.0]
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
     
     absTol=1e-5
     relTol=1e-2   
     sol=QSS_Solve(odeprob,slvr,dQmin=absTol,detection=Val(cycleDetection),dQrel=relTol,finalTime=50.0,maxErr=10*relTol)
    # save_Sol(sol)
     solInterp=solInterpolated(sol,0.01)
     er1=getError(solInterp,1,x1)  
     er2=getError(solInterp,2,x2) 
   #  timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss1(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0,maxErr=10*$relTol)
     resmliqss1= ("$(sol.algName)",relTol,(er1+er2)/2,sol.totalSteps,sol.simulStepCount,timenmliqss)
     @show resmliqss1

     return resmliqss1

   
end
 

note="analytic_abs-5"
resmliqss10=test53(note,mliqss1_SimulAna(),0)
resmliqss11=test53(note,mliqss1_SimulAna(),1)
resmliqss12=test53(note,mliqss1_SimulAna(),2)
resmliqss13=test53(note,mliqss1_SimulAna(),3)
resmliqss14=test53(note,mliqss1_SimulAna(),4)

resmliqss20=test53(note,mliqss2_SimulAna(),0)
resmliqss21=test53(note,mliqss2_SimulAna(),1)
resmliqss22=test53(note,mliqss2_SimulAna(),2)
resmliqss23=test53(note,mliqss2_SimulAna(),3)
resmliqss24=test53(note,mliqss2_SimulAna(),4)

resmliqss30=test53(note,mliqss3_SimulAna(),0)
resmliqss31=test53(note,mliqss3_SimulAna(),1)
resmliqss32=test53(note,mliqss3_SimulAna(),2)
resmliqss33=test53(note,mliqss3_SimulAna(),3)
resmliqss34=test53(note,mliqss3_SimulAna(),4)

XLSX.openxlsx("sys53_solvers_$note.xlsx", mode="w") do xf
  sheet = xf[1]
  sheet["A1"] = "sys53_solvers__$note)"
  sheet["A4"] = collect(("solver","Tolerance","error","totalSteps","simul_steps","time"))
  sheet["A5"] = collect(resmliqss10)
  sheet["A6"] = collect(resmliqss11)
  sheet["A7"] = collect(resmliqss12)
  sheet["A8"] = collect(resmliqss13)
  sheet["A9"] = collect(resmliqss14)

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