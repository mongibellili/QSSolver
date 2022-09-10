using qss
using BenchmarkTools

using Plots;
gr();
function test()
    Debug = false
    psettings = ProblemSettings(5.0, saveat(0.1), qss2())
    odeprob = @NLodeProblem begin
        u = [1.0, 2.4]
        d = [1.0, 0.5]
        dq[1] = q[2]#+2.0q[2]*q[1]   
        dq[2] = -q[1] - q[2]
        # du[2]=u[2]-u[1]+cos(t)
        #=    du1=u1+2.0u2   
           du2=u2-3#*cos(t) =#



        #=     

             du1=u3-u2*u1+2.0u2     # du1...2...N are expected to be in order....later can fix this
             du2=-u1-3u2-d1
             du3=u3+u2+d1+d2   
             if u1+0.7+d2 >0   #still have not added 'usercode checking'....throw error msg
                 d1=0.1
             else
                 d1=1.0
                 u2=4.5
             end
             if u2+d1 >0 
                 d2=0.33
                 u3=2.2
             else
                 d2=1.0
             end =#

    end

    #eqs=module1.prepareProb(odeprob)
    #@btime module1.integrate(u,eqs)
    sol = QSS_Solve(psettings, odeprob)
    #display(sol)
    if !Debug
        temp1 = []
        temp2 = []
        #if sol isa Tuple
            for i = 1:length(sol[2][1])
                push!(temp1, sol[2][1][i].coeffs[1])
                push!(temp2, sol[2][2][i].coeffs[1])
            end
            display(plot!(sol[1], temp1))
            display(plot!(sol[1], temp2))
       # end
    end



end

#@btime test()
test()
#=  using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2]#*cos(t) #*30*exp(t)
        du[2] = -u[1] - u[2]#+sqrt(t) # +1#+10-t*t+t +cos(t)
       # du[2]=(1 - u[1]^2) * u[2] - u[1] 
      # du[2]=1/(t+1) + sin(t)*sqrt(t)
     # du[2]=sqrt(t)

    end
    tspan = (0.0,5)
    u0 = [1.0,2.4]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end
odeDiffEquPackage() =#



#@btime modqss.qss_Integrate(initCond,jacobian,order)


println("done")
readline()
