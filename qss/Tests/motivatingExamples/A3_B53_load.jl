

using formalqss

#include("/home/unknown/formalqss/models/BProblem.jl")
#= function B53(j::Int, zc::Int, ev::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        subadd(mulT(-20.0, q[1], cache[2]), mulT(80.0, q[2], cache[3]), 1600.0, cache[1])
        return nothing
    elseif j == 2
        subadd(mulT(1.24, q[1], cache[2]), mulT(0.01, q[2], cache[3]), 0.2, cache[1])
        return nothing
    end
end =#
include("/home/unknown/formalqss/Tests/trash/testSave2.jl")

function test53()
    #=  pr= @saveNLodeProblem begin
         sysb53;
         "/home/unknown/formalqss/Tests/trash/testSave.jl";
         u = [-1.0, -2.0]
        # discrete = [0.0]
         du[1] = -20.0*u[1]-80.0*u[2]+1600.0
         du[2] =1.24*u[1]-0.01*u[2]+0.2
     end  
     @show pr =#
     
     u1, u2 = -8.73522174738572, -7.385745994549763
     λ1, λ2 = -10.841674966758294, -9.168325033241706
     c1, c2 = 121.14809142478035, -143.14809142478035
     xp1, xp2 = 0.0, 20.0
     x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
     x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
     solmliqss2=QSS_Solve_from_model(prob_advection10d001(),100.0,qss2(),saveat(0.01),0.0,1e-6,1e-3)
     solmliqss2Interp=solInterpolated(solmliqss2,0.1,100.0)
     #save_Sol(solmliqss2,"x1",1;xlims=(30.0,40.0),ylims=(-0.000005,0.000005))
     er1=getError(solmliqss2Interp,1,x1)  
     er2=getError(solmliqss2Interp,2,x2) 
     resnmliqss= ("$(solmliqss2.algName)",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
     @show resnmliqss

  end

  test53()
  
  