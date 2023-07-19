using OrdinaryDiffEq
#using BenchmarkTools
using BSON
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
      d=1.0
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

    solFeagin14 = solve(prob,Feagin14(),saveat=0.01,abstol = 1e-12, reltol = 1e-8)

     solFeagin14VectorN1000d1=solFeagin14.u

     BSON.@save "formalqss/ref_bson/solVectAdvection_N1000d1_Feagin14e-12.bson" solFeagin14VectorN1000d1

#= println("breakpt1")

      BSON.@load "formalqss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
      
      println("breakpt2")  =#


    #=   solTsit5= solve(prob,Tsit5(),saveat=0.01,abstol = 1e-9, reltol = 1e-5)
      println("breakpt3") 
      err2=getAverageErrorByRefs(solTsit5.u,solFeagin14VectorN1000d01,10,1000)
 @show err2
 @btime solve($prob,Tsit5(),saveat=0.01,abstol = 1e-9, reltol = 1e-5)  =#

#= solRodas5P = solve(prob,BS3(),saveat=0.01,abstol = 1e-9, reltol = 1e-5) #1.235 ms (1598 allocations: 235.42 KiB)
println("breakpt3") 
err2=getAverageErrorByRefs(solRodas5P.u,solFeagin14VectorN1000d01,10,1000)
 @show err2

 @btime solve($prob,BS3(),saveat=0.01,abstol = 1e-9, reltol = 1e-5)  =#

#=  println("breakpt3")
 solQNDF2 = solve(prob,QNDF2(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) #
 println("breakpt4")
 err1=getAverageErrorByRefs(solQNDF2.u,solFeagin14VectorN1000d01,10,1000) =#

#=  @show err1
 @btime solve($prob,QNDF2(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) 
 println("breakpt5") =#



#=  solRosenbrock23 = solve(prob,Rosenbrock23(),saveat=0.01,abstol = 1e-10, reltol = 1e-6) #
 println("breakpt3") 
 err2=getAverageErrorByRefs(solRosenbrock23.u,solFeagin14VectorN1000d01,10,1000)
 @show err2 =#
 #= solQBDF2 = solve(prob,QBDF2(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) #
 err3=getAverageErrorByRefs(solQBDF2.u,solFeagin14VectorN1000d01,10,1000)
 @show err3 =#
 #@btime solve($prob,TRBDF2(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) 
 #@btime solve($prob,Rosenbrock23(),saveat=0.01,abstol = 1e-10, reltol = 1e-6) 
 #@btime solve($prob,QBDF2(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) 



end
#@time 
odeDiffEquPackage()  


