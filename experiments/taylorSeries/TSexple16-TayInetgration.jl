# not to be called along "using OrdinaryDiffEq"
using Plots;gr()
using BenchmarkTools

#function taylorodeDiffEquPackage()

using TaylorIntegration 
     function henonheiles!(dx, x, p, t)
        dx[1] = x[2]
       # dx[2] = -x[1]-x[2]
        dx[2] = 1/(t+1) + sin(t)*sqrt(t)
        nothing
    end
    tv = 0.0:0.1:5.0
    #tv2=(0.0,5)
    x0=[1.0,2.0]
       xv = taylorinteg(henonheiles!, x0, tv, 3, 1e-6, maxsteps=30000 )
       s1=Int(length(xv)/2)
      # display(s1)
   #display(xv[1:s1])
   #display(xv[s1+1:end])
   display(plot!(tv, xv[1:s1]))
   display(plot!(tv, xv[s1+1:end])) 
#end
#= tT, xT = taylorinteg(diffeq, 3.0, 0.0, 0.34, 28, 1e-7, maxsteps=166);



title!("x vs t (log-log)")
xlabel!(L"\log_{10}(t)")
ylabel!(L"\log_{10}(x(t))") =#

#@time taylorodeDiffEquPackage()

using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2] #+t+1
        #du[2] = -u[1] - abs(u[2]) # +1#+10-t*t+t +cos(t)   
        du[2]=1/(t+1) + sin(t)*abs(2-t)
    end
    tspan = (0.0,5)
    u0 = [1.0,2.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
    #display(sol)
    display(plot!(sol,line=(:dot, 4)))
end
 #odeDiffEquPackage()


println("done") 
readline()
