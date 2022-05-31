using TaylorSeries
displayBigO(false)


println("------------Tutorial 4--------------")
#problem: solve this ODE dx/dt=x*x=f
f(x) = x*x
#the solution is going to be a taylor expansion of the exacte sol
degree = 20
u0= Taylor1([1],degree)#u0 is our initial guessed expansion. here we guessed 1.0
∫⬩dt(u::Taylor1) = integrate(u)
function taylor_step(f, u0)
    u = copy(u0)
    println("u: ",u)
    unew = u0 + ∫⬩dt(f(u))
   #  println("unew: ",unew)
    while unew != u
        u = unew

        unew = u0 + ∫⬩dt(f(u))   # Picard iteration
     #   println("unew: ",unew)
    end 
    return u
end

picardTaylorSol=taylor_step(f, u0)
@show picardTaylorSol

println("--compare with exacte sol--")

g(a) = 1/(1-a)
#@show g(u1) for taylor1([1], 3) error: division does not define a talor poly
#taylor expansion is not allowed for franctions
#u1= Taylor1([2],degree)
#exacteTayloredSol=g(u1)
#@show exacteTayloredSol
#@show exacteTayloredSol(2.3)
#the following displays show that the expansion only valid between -1 and 1
@show picardTaylorSol(0.1)
@show g(0.1)
@show picardTaylorSol(0.9)
@show g(0.9)
@show picardTaylorSol(-0.1)
@show g(-0.1)
@show picardTaylorSol(-0.5)
@show g(-0.5)
@show picardTaylorSol(1.1)
@show g(1.1)
