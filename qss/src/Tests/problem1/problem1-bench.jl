
using BenchmarkTools
using qss
function test()
    odeprob = @NLodeProblem begin
        cacheT=10 #recommend the cache to match the largest equation you have using the specs of "caches/ops"
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]
    end
   sol = QSS_Solve(odeprob,5.0,qss2())

end
@btime test()
