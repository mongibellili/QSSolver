using BenchmarkTools
using DifferentialEquations
#using Plots
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
BSON.@load "test/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
  
function test()
    function adr(du,u,p,t)# api requires four args
        _dx=100.0
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
    #Construct the problem
    odeprob = ODEProblem(adr,u0,tspan)
    #Solve the problem
    sol=solve(odeprob,ImplicitEuler()#= ,saveat=0.01 =#,abstol=1e-2,reltol=1e-1)
    #@show sol.stats
    @show length(sol.u)
   err2=getAverageErrorByRefs(sol.u,solFeagin14VectorN1000d01,1000,1000)
    @show  err2
    @btime solve($odeprob,ImplicitEuler(),abstol=1e-3,reltol=1e-2)
    #save_Sol(sol,1,400,1000)
 #=    p1=plot(sol,idxs=[1,400,1000],title="ImplicitEuler");
    savefig(p1, "adr_ImplicitEuler")  =#
    end
    #@btime 
    test()