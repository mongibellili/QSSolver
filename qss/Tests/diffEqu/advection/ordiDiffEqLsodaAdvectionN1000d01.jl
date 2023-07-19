using OrdinaryDiffEq
using LSODA
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

println("breakpt1")

      BSON.@load "formalqss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
      
      println("breakpt2") 

  aTol=1e-5
  relTol=1e-2
  sollsoda = solve(prob,lsoda(),saveat=0.01,abstol = aTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
  println("breakpt4") 
  err2=getAverageErrorByRefs(sollsoda.u,solFeagin14VectorN1000d01,1000,1000)
  @show  err2
  ttlsoda=@belapsed solve($prob,lsoda(),saveat=0.01,abstol = $aTol, reltol = $relTol) 
  reslsoda1= ("lsoda",err2,ttlsoda)
  @show reslsoda1
  

  aTol=1e-6
  relTol=1e-3
  sollsoda = solve(prob,lsoda(),saveat=0.01,abstol = aTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
  println("breakpt4") 
  err2=getAverageErrorByRefs(sollsoda.u,solFeagin14VectorN1000d01,1000,1000)
  @show  err2
  ttlsoda=@belapsed solve($prob,lsoda(),saveat=0.01,abstol = $aTol, reltol = $relTol) 
  reslsoda2= ("lsoda",err2,ttlsoda)
  @show reslsoda2

  aTol=1e-7
  relTol=1e-4
  sollsoda = solve(prob,lsoda(),saveat=0.01,abstol = aTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
  println("breakpt4") 
  err2=getAverageErrorByRefs(sollsoda.u,solFeagin14VectorN1000d01,1000,1000)
  @show  err2
  ttlsoda=@belapsed solve($prob,lsoda(),saveat=0.01,abstol = $aTol, reltol = $relTol) 
  reslsoda3= ("lsoda",err2,ttlsoda)
  @show reslsoda3




 XLSX.openxlsx("ADR_diffEq_lsoda_N1000d01_.xlsx", mode="w") do xf
  sheet = xf[1]
  sheet["A1"] = "ADR_diffEq_implicit_N1000d01 lsoda__"
  sheet["A2"] = "__e-2_e-5)"
  sheet["A3"] = collect(("solver","error","time"))
  sheet["A4"] = collect(reslsoda1)
  sheet["A5"] = "__e-3_e-6)"
  sheet["A6"] = collect(("solver","error","time"))
  sheet["A7"] = collect(reslsoda2)
  sheet["A8"] = "__e-4_e-7)"
  sheet["A9"] = collect(("solver","error","time"))
  sheet["A10"] = collect(reslsoda3)
  #= for k=1:1000
    sheet["A$(3*k+4)"]="variable $k changes: $(length(solnmliqss2.savedTimes[k]))"
    sheet["A$(3*k+5)"]=collect(solnmliqss2.savedTimes[k])
    sheet["A$(3*k+6)"]=collect(solnmliqss2.savedVars[k])
  end 
  =#
end

end
#@time 
odeDiffEquPackage()  


