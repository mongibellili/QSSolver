

using formalqss

using BenchmarkTools

include("D:/models/BProblem.jl")

function test53(#= ft::Float64,solver::Val{V} =#)
    
    odeprob = @NLodeProblem begin
     #sys b53
    u = [-1.0, -2.0]
    discrete = [0.0]
    du[1] = -20.0*u[1]-80.0*u[2]+1600.0
    du[2] =1.24*u[1]-0.01*u[2]+0.2
end  
@btime solmliqss2=QSS_Solve_from_model(B53,$odeprob,100.0,liqss2(),saveat(0.01),0.0,1e-6,1e-3)
@btime solmliqss2=QSS_Solve_from_model(B53,$odeprob,100.0,mliqss2(),saveat(0.01),0.0,1e-6,1e-3)
@btime solmliqss2=QSS_Solve_from_model(B53,$odeprob,100.0,nliqss2(),saveat(0.01),0.0,1e-6,1e-3)
@btime solmliqss2=QSS_Solve_from_model(B53,$odeprob,100.0,nmliqss2(),saveat(0.01),0.0,1e-6,1e-3)
    
     
  end

  test53()
  