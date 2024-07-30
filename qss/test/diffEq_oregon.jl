

using OrdinaryDiffEq
#= using BenchmarkTools
using XLSX =#
using Plots



function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
        du[1] = 100.8*(9.523809523809524e-5*u[2]-u[1]*u[2]+u[1]*(1.0-u[1]))
        du[2] =40320.0*(-9.523809523809524e-5*u[2]-u[1]*u[2]+u[3])
        du[3] = u[1] -u[3]
    end
    tspan = (0.0,25.0)
    u0= [1.0, 1.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)


   absTol=1e-3
   relTol=1e-2

   solRodas5P = solve(prob,Rodas5P(),saveat=0.01,abstol = 1e-12, reltol = 1e-8)
# sol = solve(prob,Rosenbrock23(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
#=  p1=plot!(sol,idxs=(0,1,2))
 savefig(p1, "plot_solRosenbrock23_vars0_1_2.png") =#

 p1=plot!(solRodas5P,marker=:circle)
 savefig(p1, "plot_sol_oregon.png")
 @show solRodas5P.stats
 solRodas5PVectorOregon=solRodas5P.u

   
 
BSON.@save "qss/ref_bson/solVectAdvection_Oregon_Rodas5Pe-12.bson" solRodas5PVectorOregon
end

odeDiffEquPackage()  


