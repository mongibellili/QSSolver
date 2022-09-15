using qss
using Plots;
gr();
include("/home/unknown/QS_Solver/qss/src/models/classicProblem.jl") # where you saved the model
function test()
    odeprob = @NLodeProblem begin
        cacheT=4 #recommend the cache to match the largest equation you have using the specs of "caches/ops"
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]
    end
  
    sol=QSS_Solve_from_model(f,odeprob,5.0,qss2()) # later if to change f to meaningfull name and maybe save all models in one file, user has to enter model name in the macro so that f be name upfront
    temp1 = []
        temp2 = []
        for i = 1:length(sol[2][1])
            push!(temp1, sol[2][1][i].coeffs[1])
            push!(temp2, sol[2][2][i].coeffs[1])
        end
        display(plot!(sol[1], temp1))
        display(plot!(sol[1], temp2))
       println("done")
        readline() 
end
test()
