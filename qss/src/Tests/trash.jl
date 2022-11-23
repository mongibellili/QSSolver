#= using Plots
#analytic to this
#= u = [1.0, 0.0]
discrete = [0.0]
du[1] = -4.0*u[1]+0.5*u[2]
du[2] =u[1]-0.25*u[2]  =#
u1=(-sqrt(257)-15)/8
u2=(sqrt(257)-15)/8
λ1=(-sqrt(257)-17)/8
λ2=(sqrt(257)-17)/8
c1=-4/sqrt(257)
c2=4/sqrt(257)
x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)
x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)
display(plot!(x1,xlims=(60,80),ylims=(0,0.0004)))
display(plot!(x2,xlims=(60,80),ylims=(0,0.0004)))
println("press enter to exit")
readline()
 =#
 savedVars=Vector{Float64}(undef, 2)
 display(savedVars);println()
 savedVars[1]=3.14
 resize!(savedVars,4)
 savedVars[3]=3.14
 display(savedVars)