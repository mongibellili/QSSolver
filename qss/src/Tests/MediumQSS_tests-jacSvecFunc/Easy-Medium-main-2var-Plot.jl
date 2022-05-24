
using StaticArrays
using qss
using Plots;gr()
 #using OrdinaryDiffEq
#using DifferentialEquations =#
function MediumqssApproachInitialInside(jacobian::SVector{2,Function}) 
    initConditions=@SVector[1.0,2.0]
   SD=@SMatrix[2 0;1 2] 
   psettings = ProblemSettings(5.0,saveat(0.1),qss1())
   prob = QSS_Problem(initConditions,jacobian,SD)  
   sol=QSS_Solve(psettings,prob)   
   #display(sol)
   display(plot!(sol[1],sol[2][1],label = ["mediumQSS sol1"]))
    display(plot!(sol[1],sol[2][2],label = ["mediumQSS sol2"]) )
end

f1(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int) where {R,T}=q[order+1]  
f2(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int) where {R,T}= -q[1]-q[order+1]
jacobian=SVector{2,Function}(f1,f2)  
#@btime MediumqssApproachInitialInside(jacobian) 
MediumqssApproachInitialInside(jacobian) 




using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2] #+t+1
        du[2] = -u[1] - u[2]  #+10-t*t+t +cos(t)   
    end
    tspan = (0.0,5)
    u0 = [1.0,2.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4),label = ["diffequ sol1" "diffequ sol2"]))
end
odeDiffEquPackage()




println("done") 
readline()









#A=[3 -2; 2 -2]; initcond=[1,1]
#anal sol1=(2/3) *exp(2x)+(1/3) *exp(-x)
#analy sol2=(1/3) *exp(2x)+(2/3) *exp(-x)
#= f(x) = (1/3) *exp(2x)+(2/3) *exp(-x)
display(plot!(f, 0, 3))=#
