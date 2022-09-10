using TaylorSeries
#use_show_default(true)
#= q1=Taylor1([1.1,1.0],4)
println(q1)
#= v1=Taylor1([2.0, 0.0,0.0])
println(q1(v1))
v11=Taylor1([0.0, 2.0,1.0])
println(q1(v11))
v111=Taylor1([2.0, 2.0])
println(q1(v111)) =#
v1111=Taylor1([2.0, 1.0])
println(q1(v1111))
println(q1) =#
#= q1=Taylor1([2.0,1.3,1.0],4)
v1=Taylor1([1.1, 0.0])
println(q1(v1))
q1=Taylor1([2.0,1.3,1.0],4)
v1=Taylor1([1.1, 0.0])
println(q1(v1))
 =#
t=3.0
 q=Taylor1([2.0,3.0,1.0],3)
f(q,t)=q
g(q,t)=q*cos(t)
#println(differentiate(f(q,t))*cos(t)==differentiate(g(q,t)))
t2=Taylor1(3)
q=Taylor1([2.0,3.0,1.0],3)
#println(differentiate(f(q,t2))*cos(t2)==differentiate(g(q,t2)))


println((differentiate(g(q,t))))
println((differentiate(g(q,t2))))
println((differentiate(g(q,t)))(t))
println((differentiate(g(q,t2)))(t))   # as expected cos(t) gets evaluated right away before it gets differentiated which defeats the whole purpose of d(dx/dt)dt==> variable t has to be type taylor1