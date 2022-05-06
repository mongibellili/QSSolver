using BenchmarkTools
using StaticArrays
using qss
function qssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0,1.0]
    jacobian=@SMatrix[0.0 1.0 -1.0;-1.0 -1.0 0.0;1.0 0.0 -2.0]
    settings = ModelSettings(initConditions,jacobian,5.0,saveat(0.2),qss1())#do not call saveat to not save, i should fix when called with zero; it does not save at all
    sol=QSS_simGenerate(settings)
end
@btime qssApproachInitialInside()






