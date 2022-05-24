 using BenchmarkTools
using StaticArrays
using qss
 #using OrdinaryDiffEq
#using DifferentialEquations =#
function MediumqssApproachInitialInside(jacobian::SVector{2,Function}) 
    initConditions=@SVector[1.0,2.0]
   SD=@SMatrix[2 0;1 2] 
   psettings = ProblemSettings(1.0,saveat(0.1),qss1())
   prob = QSS_Problem(initConditions,jacobian,SD)  
   sol=QSS_Solve(psettings,prob)   
end



#= order=1
R=2
T=2s
f1(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int)=q[order+1]  
f2(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int) = -q[1]-q[order+1] =#
f1(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int) where {R,T}=q[2]  
f2(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int) where {R,T}= -q[1]-q[2]+1

jacobian=SVector{2,Function}(f1,f2)  
println(typeof(jacobian[1]))
@btime MediumqssApproachInitialInside(jacobian) 
#MediumqssApproachInitialInside(jacobian) 











#A=[3 -2; 2 -2]; initcond=[1,1]
#anal sol1=(2/3) *exp(2x)+(1/3) *exp(-x)
#analy sol2=(1/3) *exp(2x)+(2/3) *exp(-x)
#= f(x) = (1/3) *exp(2x)+(2/3) *exp(-x)
display(plot!(f, 0, 3))=#
