 using BenchmarkTools
using StaticArrays
using qss

function MediumqssApproachInitialInside() 
    initConditions=@SVector[1.0,2.0]
   SD=@SMatrix[2 0;1 2] 
   psettings = ProblemSettings(5.0,saveat(0.1),qss1())
   prob = QSS_Problem(initConditions,SD)  
   sol=QSS_Solve(psettings,prob)   
end
myExp=@equations begin 
    du[1]=u[2] 
    du[2]=-u[1]-u[2] 
   # du[3]=u[3]+u[2]
end 
#= function funcX1(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int)where{R,T}
    q[2]
  end
  function funcX2(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int)where{R,T} 
   -q[1]-q[2]
  end =#

@btime MediumqssApproachInitialInside() 
#MediumqssApproachInitialInside() 











#A=[3 -2; 2 -2]; initcond=[1,1]
#anal sol1=(2/3) *exp(2x)+(1/3) *exp(-x)
#analy sol2=(1/3) *exp(2x)+(2/3) *exp(-x)
#= f(x) = (1/3) *exp(2x)+(2/3) *exp(-x)
display(plot!(f, 0, 3))=#
