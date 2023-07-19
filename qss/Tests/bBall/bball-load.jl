

using formalqss
include("/home/unknown/formalqss/Tests/trash/testSave.jl")
function test()
    
 
    sol=QSS_Solve_from_model(bball,bball_prob(),20.0,qss2(),saveat(0.01),0.0,1e-6,1e-3)
    save_Sol(sol)
 
end
#@time 
test()

