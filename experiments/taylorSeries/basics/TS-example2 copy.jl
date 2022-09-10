using TaylorSeries
using MacroTools: postwalk
using BenchmarkTools
function getcoeffsRhs(expr::Expr)
    smbl=Symbol(expr,:[r])
    str=String(smbl)
    myexpr=Meta.parse(str)
    return myexpr
end
function getcoeffsLhs(expr::Expr)
    smbl=Symbol(expr,:[l])
    str=String(smbl)
    myexpr=Meta.parse(str)
    return myexpr
end
macro tayArth2(expr)
    lhs=postwalk(x -> x isa Expr && x.head == :ref && x.args[1] != :d ? getcoeffsLhs(x) : x, expr.args[1])# if !=d and head not needed 
    #lhs=getcoeffsLhs(lhs) 
    rhs=postwalk(x -> x isa Expr && x.head == :ref && x.args[1] != :d ? getcoeffsRhs(x) : x, expr.args[2])
    computedExpr=quote
                    l=1
                    r=0
                    $lhs=$rhs
                    l=2
                    r=1
                    $lhs=$rhs/2
                           #  $lhsrhs
                            # println($lhsrhs)
                     
                end
    esc(computedExpr)
end
d=[1.2]
x=[Taylor1([0.1,0.0,0.0]),Taylor1([1.0,0.0,0.0])]
q=[Taylor1([1.1,2.2,3.0]),Taylor1([1.2,3.4,4.0])]
#= @tayArth2 x[1]=(-9.8 - (d[1] / 1.0) * (100000.0 * q[1] + 30.0 * q[2])) =#

function test(d::Vector{Float64},x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}})
    @tayArth2 x[1]=(-2.0 - (d[1] / 1.0) * (3.0 * q[1] + 2.0 * q[2]))
end
#@btime 
test(d,x,q)
@show x