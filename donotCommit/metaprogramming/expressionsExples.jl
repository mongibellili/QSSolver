#=function math_expr(op, op1, op2)
    #expr=Expr(:call, op, op1,op2)
    #expr=:(:call, op, op1,op2)
    symbl=Symbol(op1,op,op2)
    str=String(symbl)
    expr=Meta.parse(str)
    return (expr)
      
end
ex = math_expr(+, 1, 5)
@show eval(ex)
=#
#import Base: *
function math_expr2(op, op1, op2)
    (op1N,op2N)=map(x -> isa(x,Number) ? 2*x : x,(op1,op2))
    expr=Expr(:call, op, op1N, op2N)
    return expr
      
end

Base.:*(x::String, y::Int64) = string(x, y) # base is needed here, unlike in promotion.jl cuz eval expr goes for Base.:*
#@show "3"*2
ex2 = math_expr2(*, "e", 2)
@show eval(ex2)