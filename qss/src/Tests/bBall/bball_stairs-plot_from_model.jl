

using qss
using Plots;
include("/home/unknown/QS_Solver/qss/src/models/BBall.jl")# the location where you saved
function test()
    odeprob = @NLodeProblem begin
        k=1e+5
        g=9.8
        b=30.0
        u = [10.5,0.0,0.575,0.5]
        discrete = [0.0,10.0]

        #y=1 vy=2 x=3 vx=4
        du[1] = u[2];
        d[2] = -9.8 - 0.1 * u[2] - discrete[1] * ((u[1] - discrete[2]) *1e6+ u[2] * 30);
        d[3] = u[4];
        d[4] = -0.1 * u[4];

        du[1] =u[2]
        du[2] =-g-(discrete[1])*(k*u[1]+b*u[2])
        if -u[1]+discrete[2]>0   
            discrete[1]=1.0               
        else
            discrete[1]=1.0                                    
        end
        if u[3]-11+discrete[2]>0
            discrete[2]=discrete[2]-1
        else
        end
    end
    sol=QSS_Solve(odeprob,8.0,qss2())
    x=evaluateSol(sol,1,1.0)
    pointt=[1.0]
    pointx=[x]
    display(scatter(pointt, pointx,label="x1",line=(:dot, 6),color=:black))
    plotSol(sol)
end
test()
