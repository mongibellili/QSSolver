
#using Plots;gr()
using BenchmarkTools
using StaticArrays
using qss
#using TimerOutputs  
function qssApproach(initConditions::SVector{N,Float64},jacobian::SMatrix{2,2,Float64,4}) where {N}
#=     initialTime=0.0
    finalTime=0.5 # right now it only accepts float
    dQmin=1e-6
    dQrel=1e-3
    order=1 # order2 means we will consider second derivatives =#
    #solver=qss1()

   # initConditions=@SVector[1.0,2.0]
 
   # jacobian=@SMatrix[0.0 1.0;-1.0 -1.0 ]
   # display(jacobian)

   settings = ModelSettings(initConditions,jacobian,0.5,0.0,1e-6,1e-3,1,Val(1),0.1)
   # settings = ModelSettings(initConditions,jacobian,finalTime,saveat(0.1))
    #display(typeof(settings.initConditions))
   # QSS_integrate(settings)
    #display(simulator)
    #plotX(simulator)
end
#initConditions=[1.0,2.0]
 initConditions=@SVector[1.0,2.0]
 #jacobian=[0.0 1.0;-1.0 -1.0 ]
 jacobian=@SMatrix[0.0 1.0;-1.0 -1.0 ]
#qssApproach(initConditions,jacobian) 
 @btime qssApproach(initConditions,jacobian)
