

using formalqss
#using BenchmarkTools
#= using Plots;
gr(); =#
#include("D:/models/bball.jl")
function test()
    odeprob = @saveNLodeProblem begin
        bball;
        "/home/unknown/formalqss/Tests/trash/testSave.jl";
        #= parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001 =#
        u = [20.0,0.0]
        discrete = [15.0]
        du[1] =u[2]
        du[2] =-9.8
        if discrete[1]-u[1]>0   #5*discrte gave error
            u[2]=-u[2]   #discrete=0.0-->type Symbol has no field args...find to personalize error msg  
                   
       #=  else
            discrete[1]=-1.0   =#                                  
        end
        if u[2]>0
            discrete[1]=discrete[1]-5.0
        end
     #=    if u[2]>0
            u[1]=5.0

        end =#
    end
   # @show odeprob
 
   #=  sol=QSS_Solve_from_model(bball,odeprob,10.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(sol) =#
 
end
#@time 
test()
#= function test2()
    odeprob = @NLodeProblem begin
        parameter2=30.0# why 300 bad
        parameter1=1e5
        u = [10.0,0.0]
        discrete = [0.0]
        du[1] =u[2]
        du[2] =-9.8-(discrete[1])*(parameter1*u[1]+parameter2*u[2])
        if -u[1]>0   
            discrete[1]=1.0           
        else
            discrete[1]=0.0                                    
        end
    end
    sol= QSS_Solve(odeprob,1.5,qss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(sol)
    
   #=  sol= QSS_Solve(odeprob,1.0,qss3(),saveat(0.01),0.0,1e-6,1e-3)
      save_Sol(sol) =#
end
test2() =#
