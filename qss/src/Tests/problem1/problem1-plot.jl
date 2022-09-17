using qss
using Plots;
function test()
    odeprob = @NLodeProblem begin
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]
    end
   sol = QSS_Solve(odeprob,2.0,qss2())
#=   x=evaluateSol(sol,1,1.0)
  pointt=[1.0]
  pointx=[x]
  display(scatter(pointt, pointx,label="x1",line=(:dot, 6),color=:black))
   plotSol(sol) =#


   
end
test()
