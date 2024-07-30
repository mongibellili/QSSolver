using OrdinaryDiffEq
#= using BenchmarkTools
using XLSX =#
using Plots



function odeDiffEquPackage()
 
    function funcName(du,u,p,t)# api requires four args
       #=  du[1] = acos(sin(u[2]))
        du[2] = (u[1]) =#
        #= du[1] = t
        du[2] =1.24*u[1]-0.01*u[2]+0.2 =#
       
       #=  du[1] = t
        
        for k in 2:5 
            du[k]=(u[k]-u[k-1]) ;
        end   =#
        du[1] = -20.0*u[1]-80.0*u[2]+1600.0
         du[2] =1.24*u[1]-0.01*u[2]+0.2
     #=    du[1] = -u[2]+u[1]
    du[2]=u[1]-u[2] =#
    end
    tspan = (0.0,0.5)
    u0= [-1.0, -2.0]
    prob = ODEProblem(funcName,u0,tspan)


   absTol=1e-3
   relTol=1e-2

sol = solve(prob,ImplicitEuler(),abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
# sol = solve(prob,Rosenbrock23(),saveat=0.01,abstol = absTol, reltol = relTol) #1.235 ms (1598 allocations: 235.42 KiB)
#=  p1=plot!(sol,idxs=(0,1,2))
 savefig(p1, "plot_solRosenbrock23_vars0_1_2.png") =#
 @show length(sol.t)
 p1=plot!(sol,marker=:circle)
 savefig(p1, "plot_sol_ft05.png")
 @show sol.stats
end

odeDiffEquPackage()  



