using TaylorSeries
using BenchmarkTools
displayBigO(false)


println("------------Tutorial 5--------------")
#problem: solve this ODE dx/dt=x*x=f where x is vector and f R^d --> R^d 

f(x::Vector) = x.*x
#the solution is going to be a taylor expansion of the exacte sol
#v=[2,3,4]
#display(v)
#println()
#display(f(v))
degree = 3
u0= Taylor1([1],degree) #u0 is our initial guessed expansion. here we guessed 1.0
v0=[u0,u0,u0]
#display(v)
#println()
#display(f(v))
∫⬩dt(u::Taylor1) = integrate(u)
function taylor_step(f, v0)
    u = copy(v0)
    unew = copy(v0)
    #println("u: ",u)
    #map!((a,b)->a+∫⬩dt.(f(b)),unew,v0,u)
    unew = v0 .+ ∫⬩dt.(f(u))
    #println("unew: ",unew)
    while unew != u
        u = unew

        unew = v0 .+ ∫⬩dt.(f(u))   # Picard iteration
       # println("unew: ",unew)
    end 
    return u
end

#@show taylor_step(f, v0)

@btime taylor_step(f, v0)


