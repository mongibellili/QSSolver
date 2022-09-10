using ApproxFun
using BenchmarkTools


#= a = Fun(cos,Chebyshev())
display(a.coefficients) =#
#= a = Fun(cos,Taylor())
#display(a.coefficients)
display(a(1.57)) =#
#= a = Fun(sqrt,Taylor()) #Maximum number of coefficients 1048577 reached in constructing Fun.
display(a) =# 

#= x = Fun(identity,-1..1);
display(exp(x));println()
display(cos(x));println()
display(sqrt(x));println() =#
function test()
    x = Fun(identity,-1..1);
    for i=1:10
        exp(x+i)
    end
end
@btime test()