using TaylorSeries
using MacroTools: postwalk
using BenchmarkTools
function getcoeffs(expr::Expr)
    smbl=Symbol(expr,:[k])
    str=String(smbl)
    myexpr=Meta.parse(str)
    return myexpr
end
macro tayArth2(expr)
    Base.remove_linenums!(expr)
    lhsrhs=postwalk(x -> x isa Expr && x.head == :ref && x.args[1] != :d ? getcoeffs(x) : x, expr)
    #println(lhsrhs)
    computedExpr=quote
                    for k=0:2
                             $lhsrhs
                            # println($lhsrhs)
                     end
                end
    esc(computedExpr)
end
d=[1.0]
x=[Taylor1([0.0,0.0,0.0]),Taylor1([0.0,0.0,0.0])]
q=[Taylor1([1.0,2.0,3.0]),Taylor1([4.0,5.0,6.0])]
#= @tayArth2 x[1]=(-9.8 - (d[1] / 1.0) * (100000.0 * q[1] + 30.0 * q[2])) =#

function fasttest(d::Vector{Float64},x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}})
   # @tayArth2 x[1]=(-2.0 +  (3.0 * q[1] + 2.0 * q[2]))
    @tayArth2 x[1] = (-2.0 +3*q[1]+   2.0 *q[2])# investigate error (::LineNumberNode, ::Module, ::Any)
    #12.556 ns (0 allocations: 0 bytes)
end
function test(d::Vector{Float64},x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}})
   # x[1]=(-2.0 +  (3.0 * q[1] + 2.0 * q[2]))
   x[1]=-2.0 +3*q[1] + 2.0 *q[2]#174.420 ns (4 allocations: 320 bytes)
end
#@btime fasttest(d,x,q) 
#fasttest(d,x,q)
#@show x
println(@macroexpand @tayArth2 x[1] = (-2.0 +3*q[1]+   2.0 *q[2]))