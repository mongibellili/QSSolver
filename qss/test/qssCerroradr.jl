using DelimitedFiles
using DifferentialEquations


function getAverageErrorByRefs(solref,savedTimes::Vector{Vector{Float64}},savedVars :: Vector{Vector{Float64}},T::Int) 
    allErrors=0.0
    for index=1:T
        numPoints=length(savedTimes[index])
        sumTrueSqr=0.0
        sumDiffSqr=0.0
        relerror=0.0
        for i = 1:numPoints
            timept=savedTimes[index][i]
            Ns=savedVars[index][i] #numericl sol
            ts=solref(timept,idxs=index) #true sol or reference sol
            sumDiffSqr+=(Ns-ts)*(Ns-ts)
            sumTrueSqr+=ts*ts
        end
        if  abs(sumTrueSqr)>1e-9
        relerror=sqrt(sumDiffSqr/sumTrueSqr)
        else
          relerror=0.0
        end
        #@show relerror,index
        allErrors+= relerror
    end
    return allErrors/T
  end

  


function test()
        T=1000
        savedTimes=Vector{Vector{Float64}}(undef, T)
        savedVars = Vector{Vector{Float64}}(undef, T)
        for k in 1:T
            savedTimes[k]=Vector{Float64}()
            savedVars[k]=Vector{Float64}()
            data1 = readdlm("qss/Csolver/mliqss_adr/u[$k].dat")
            for i in 1:size(data1, 1)
                push!(savedTimes[k],data1[i, 1])
                push!(savedVars[k],data1[i, 2])
            end
        end
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
        solref=solve(odeprob,Feagin14(),abstol = 1e-12, reltol = 1e-8)
        err=getAverageErrorByRefs(solref,savedTimes,savedVars,T) 
        @show err

end

test()