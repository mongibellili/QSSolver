using qss
using BenchmarkTools
#= using Plots;
gr(); =#

function test()
    odeprob = @NLodeProblem begin
        #= parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001 =#
        u = [20.0,0.0]
        discrete = [15.0]
        du[1] =cos(t+u[2])
        du[2] =-9.8
        if t-5.0>0.0   #5*discrte gave error
            u[2]=-u[2]  #discrete=0.0-->type Symbol has no field args...find to personalize error msg  
           # discrete[1]=discrete[1]-5.0          
       #=  else
            discrete[1]=-1.0   =#                                  
        end
    end
   sol= QSS_Solve(odeprob,qss2()#= ,dQmin=absTol,saveat=0.01,dQrel=relTol =#,finalTime=15.0)
   #save_Sol(sol)
   save_Sol(sol,xlims=(0.0,15.0) ,ylims=(19.0,22.0))
end
#@time 
test()

