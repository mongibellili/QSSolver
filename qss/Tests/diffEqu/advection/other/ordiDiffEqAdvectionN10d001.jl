using OrdinaryDiffEq
using BSON
using BenchmarkTools


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
   
      _dx=1.0#1/dx=N/10=10/10
        a=1.0
        d=0.001
        r=1000.0
       # discrete = [0.0]
        du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
        
        for k in 2:9  
            du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
        end 
        du[10]=-a*_dx*(u[10]-u[9])+d*_dx*_dx*(2.0*u[9]-2.0*u[10])+r*u[10]*u[10]*(1.0-u[10]) 
        
    

    
    end
    tspan = (0.0,10.0)
    u0= [1.0, 1.0, 1.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)


   #=  solFeagin14 = solve(prob,Feagin14(),saveat=0.01,abstol = 1e-12, reltol = 1e-8)

    solFeagin14VectorN10d001=solFeagin14.u

    @show length(solFeagin14VectorN10d001)
    @show solFeagin14VectorN10d001[1]
    @show solFeagin14VectorN10d001[end]
 
   BSON.@save "solVectAdvection_N10d001_Feagin14e-12.bson" solFeagin14VectorN10d001
 =#
 BSON.@load "/home/mongi/projects/formalqss/solVectAdvection_N10d001_Feagin14e-12.bson" solFeagin14VectorN10d001
 solFTRBDF2 = solve(prob,TRBDF2(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) #1.235 ms (1598 allocations: 235.42 KiB)
 err1=getAverageErrorByRefs(solFTRBDF2.u,solFeagin14VectorN10d001,10,1000)
 @show err1
 solRosenbrock23 = solve(prob,Rosenbrock23(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) #1.235 ms (1598 allocations: 235.42 KiB)
 err2=getAverageErrorByRefs(solRosenbrock23.u,solFeagin14VectorN10d001,10,1000)
 @show err2
 solQBDF2 = solve(prob,QBDF2(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) #1.235 ms (1598 allocations: 235.42 KiB)
 err3=getAverageErrorByRefs(solQBDF2.u,solFeagin14VectorN10d001,10,1000)
 @show err3
 @btime solve($prob,TRBDF2(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) 
 @btime solve($prob,Rosenbrock23(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) 
 @btime solve($prob,QBDF2(),saveat=0.01,abstol = 1e-5, reltol = 1e-5) 
end
#@btime 
odeDiffEquPackage()  



