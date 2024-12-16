using DifferentialEquations
using BSON
#using Plots
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
  BSON.@load "qss/ref_bson/solRodas5PVectorTyson.bson" solRodas5PVectorTyson
function odeDiffEquPackage()
     function f(du, u, p, t)
      du[1] = u[4]-1e6*u[1]+1e3*u[2]
      du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
      du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
      du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
      du[5] = 0.015-200.0*u[2]*u[5]
      du[6] =u[4]-0.6*u[6]
    end

    u0 = [0.0,0.75,0.25,0.0,0.0,0.0]
    tspan = (0.0, 25.0)
    prob = ODEProblem(f, u0, tspan)
  


    sol = solve(prob, ImplicitEuler(),saveat=0.01, reltol=1e-2,abstol=1e-3) 
    #@show length(sol.u)
    err2=getAverageErrorByRefs(sol.u,solRodas5PVectorTyson,6,2500)
    @show err2 
 
   #=  p1=plot(sol,idxs=[5]);
    savefig(p1, "tysonBEindex5")
    p1=plot(sol);
    savefig(p1, "tysonBE")  =#
   # solRodas5PVectorTyson=sol.u

   # BSON.@save "test/solRodas5PVectorTyson.bson" solRodas5PVectorTyson 
 end

 #@btime 
 odeDiffEquPackage() 




 # to get a ref solution:

 

#= using BSON
using DifferentialEquations
#using Plots
function odeDiffEquPackage()
     function f(du, u, p, t)
      du[1] = u[4]-1e6*u[1]+1e3*u[2]
      du[2] =-200.0*u[2]*u[5]+1e6*u[1]-1e3*u[2]
      du[3] = 200.0*u[2]*u[5]-u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)
      du[4] =u[3]*(0.018+180.0*(u[4]/(u[1]+u[2]+u[3]+u[4]))^2)-u[4]
      du[5] = 0.015-200.0*u[2]*u[5]
      du[6] =u[4]-0.6*u[6]
    end

    u0 = [0.0,0.75,0.25,0.0,0.0,0.0]
    tspan = (0.0, 25.0)
    prob = ODEProblem(f, u0, tspan)
   # sol = solve(prob, Rodas5P(), reltol=1e-5,abstol=1e-6) 
   # p1=plot!(sol);
   # savefig(p1, "tyson6_") 


    sol = solve(prob, Rodas5P(), reltol=1e-8,abstol=1e-12) 
   # p1=plot!(sol);
   # savefig(p1, "tyson12_") 

    solRodas5PVectorTyson=sol.u

    BSON.@save "qss/ref_bson/solRodas5PVectorTyson.bson" solRodas5PVectorTyson 
 end
 odeDiffEquPackage()  =#