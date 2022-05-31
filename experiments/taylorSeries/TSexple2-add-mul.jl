using TaylorSeries

use_show_default(true)
#= q=Taylor1([2.0,0.0,1.0],4)
#println(q) #2.0 + 1.0 t² + 𝒪(t⁵) #Taylor1{Float64}([2.0, 0.0, 1.0, 0.0, 0.0], 4)
v=Taylor1([1.1, 0.5])
#println(v)#1.1 + 0.5 t + 𝒪(t²)  #Taylor1{Float64}([1.1, 0.5], 1)
u=q+v
println(u)
w=q*v
println(w) =#

t=Taylor1(2)
println(t)
q=Taylor1([0.75,0.0],2)#inital cond x=0.75
println(q)
f(q,t)=q*cos(t)
println(f(q,t))

println(differentiate(f(q,t)))