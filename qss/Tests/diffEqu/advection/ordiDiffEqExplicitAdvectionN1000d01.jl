using OrdinaryDiffEq
using BenchmarkTools
using BSON
using XLSX
function getAverageErrorByRefs(sol::Vector{Vector{Float64}},solRef::Vector{Any},T::Int,numPoints::Int)
  #numPoints=length(solmliqss.savedTimes[1])
  allErrors=0.0
  for index=1:T
      sumTrueSqr=0.0
      sumDiffSqr=0.0
      relerror=0.0
      for i = 1:numPoints #
          ts=solRef[i][index]
          Ns=sol[i][index]
          sumDiffSqr+=(Ns-ts)*(Ns-ts)
          sumTrueSqr+=ts*ts
      end
      relerror=sqrt(sumDiffSqr/sumTrueSqr)
      
      allErrors+= relerror
  end
  return allErrors/T
end
function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
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
    tspan = (0.0,10.0)
    
    u0=zeros(1000)
    u0[1:333].=1.0


    prob = ODEProblem(funcName,u0,tspan)

   #=  solFeagin14 = solve(prob,Feagin14(),saveat=0.01,abstol = 1e-12, reltol = 1e-8)

     solFeagin14VectorN1000d01=solFeagin14.u

     BSON.@save "formalqss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01 =#
  BSON.@load "formalqss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01   

  absTol=1e-5
    relTol=1e-2
solBS3 = solve(prob,BS3(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
println("breakpt3") 
err1=getAverageErrorByRefs(solBS3.u,solFeagin14VectorN1000d01,1000,1000)
@show  err1
 ttBS3=@belapsed solve($prob,BS3(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
 resBS3_1= ("bs3",err1,ttBS3)
 @show resBS3_1

 solRalston = solve(prob,Ralston(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
 println("breakpt4") 
 err2=getAverageErrorByRefs(solRalston.u,solFeagin14VectorN1000d01,1000,1000)
 @show  err2
  ttRalston=@belapsed solve($prob,Ralston(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
  resRalston_1= ("Ralston",err2,ttRalston)
  @show resRalston_1

  absTol=1e-6
    relTol=1e-3
solBS3 = solve(prob,BS3(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
println("breakpt3") 
err1=getAverageErrorByRefs(solBS3.u,solFeagin14VectorN1000d01,1000,1000)
@show  err1
 ttBS3=@belapsed solve($prob,BS3(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
 resBS3_2= ("bs3",err1,ttBS3)
 @show resBS3_2

 solRalston = solve(prob,Ralston(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
 println("breakpt4") 
 err2=getAverageErrorByRefs(solRalston.u,solFeagin14VectorN1000d01,1000,1000)
 @show  err2
  ttRalston=@belapsed solve($prob,Ralston(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
  resRalston_2= ("Ralston",err2,ttRalston)
  @show resRalston_2

  absTol=1e-7
  relTol=1e-4
solBS3 = solve(prob,BS3(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
println("breakpt3") 
err1=getAverageErrorByRefs(solBS3.u,solFeagin14VectorN1000d01,1000,1000)
@show  err1
ttBS3=@belapsed solve($prob,BS3(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
resBS3_3= ("bs3",err1,ttBS3)
@show resBS3_3

solRalston = solve(prob,Ralston(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
println("breakpt4") 
err2=getAverageErrorByRefs(solRalston.u,solFeagin14VectorN1000d01,1000,1000)
@show  err2
ttRalston=@belapsed solve($prob,Ralston(),saveat=0.01,abstol = $absTol, reltol = $relTol) 
resRalston_3= ("Ralston",err2,ttRalston)
@show resRalston_3

  



 XLSX.openxlsx("ADR_diffEq_explicit_N1000d01_.xlsx", mode="w") do xf
  sheet = xf[1]
  sheet["A1"] = "ADR_diffEq_implicit_N1000d01 2 solvers__"
  sheet["A2"] = "__e-2_e-5)"
  sheet["A3"] = collect(("solver","error","time"))
  sheet["A4"] = collect(resBS3_1)
  sheet["A5"] = collect(resRalston_1)
 
  sheet["A6"] = "__e-3_e-6)"
  sheet["A7"] = collect(("solver","error","time"))
  sheet["A8"] = collect(resBS3_2)
  sheet["A9"] = collect(resRalston_2)

  sheet["A10"] = "__e-4_e-7)"
  sheet["A11"] = collect(("solver","error","time"))
  sheet["A12"] = collect(resBS3_3)
  sheet["A13"] = collect(resRalston_3)

  #= for k=1:1000
    sheet["A$(3*k+4)"]="variable $k changes: $(length(solnmliqss2.savedTimes[k]))"
    sheet["A$(3*k+5)"]=collect(solnmliqss2.savedTimes[k])
    sheet["A$(3*k+6)"]=collect(solnmliqss2.savedVars[k])
  end =#
end

end
#@time 
odeDiffEquPackage()  


