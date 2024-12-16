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
        T=6
        savedTimes=Vector{Vector{Float64}}(undef, T)
        savedVars = Vector{Vector{Float64}}(undef, T)
        varnames = ["C2", "CP", "pM", "M", "Y", "yP"]
        for k in 1:T
            savedTimes[k]=Vector{Float64}()
            savedVars[k]=Vector{Float64}()
            data1 = readdlm("qss/Csolver/mliqss_TYSON-2-5/$(varnames[k]).dat")
            for i in 1:size(data1, 1)
                push!(savedTimes[k],data1[i, 1])
                push!(savedVars[k],data1[i, 2])
            end
        end
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

        solref = solve(prob, Rodas5P(), reltol=1e-8,abstol=1e-12) 
        err=getAverageErrorByRefs(solref,savedTimes,savedVars,T) 
        @show err

end

test()