using BenchmarkTools
using StaticArrays
using qss
function qssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0,1.0]
    inputVars=@SVector[0.0,1.0,1.0]
    jacobian=@SMatrix[0.0 1.0 -1.0;-1.0 -1.0 0.0;1.0 0.0 -2.0]
    psettings = ProblemSettings(5.0,saveat(0.1),qss1())
    prob = QSS_Problem(initConditions,jacobian,inputVars)
 
    

    sol=QSS_Solve(psettings,prob)
end
@btime qssApproachInitialInside()






