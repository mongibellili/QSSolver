using AbstractAlgebra
using BenchmarkTools


#= R, x = PowerSeriesRing(QQ, 4, "x")
@show x =#
#= f= 2 + x + 3x^2+ 3x^4
g=derivative(f)
@show f
@show derivative(f)
@show integral(g)
@show exp(x) =#
#@show cos(x)
#display(derivative(sqrt(x+0.5)))

function test()
    R, x = PowerSeriesRing(QQ, 4, "x")
    for i=1:10
        exp(x+i)
    end
end
@btime test()