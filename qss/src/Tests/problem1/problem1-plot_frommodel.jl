using qss

include("/home/unknown/QS_Solver/qss/src/models/classicProblem.jl") # where you saved the model
function test()
    odeprob = @NLodeProblem begin
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2]
    end
  
    sol=QSS_Solve_from_model(f,odeprob,5.0,qss3()) # 
    #plotSol(sol)
end
test()
