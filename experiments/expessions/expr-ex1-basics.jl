##################################exemple0#################################

#four ways to create an expression
j=1
superexpr=Meta.parse("if j==1 x=2 end")
@show superexpr
#= myex1=Meta.parse("1+5")  # from a string

#@show myex1#:(1 + 5)
 myex2=:(1+5)  #manually
 #@show myex2#:(1 + 5)
#@show myex2==myex1 #true
 #@show eval(myex)#6
 myex3=:((+)(1,5))# manually function like call
 @show myex3#:(1 + 5)
 #@show myex3==myex1 #true
 myex4=Expr(:call,+,1,5)# through the constructor. though the result is an expression that evaluates to 6 its representation is different from the others
 @show myex4#:((+)(1,5))
 # @show eval(myex4)#6
  @show myex4==myex3 #false =#
#=   myex4=Expr(:mongi)
  @show myex4
  dump(myex4)
  push!(myex4.args,1)
  @show myex4
  dump(myex4) =#

##################################exemple1#################################
                      #normal use of constructor Expr
#= function math_expr(op, op1, op2)
    expr=Expr(:call, op, op1,op2)
    @show(expr)#   :((+)(1, 5))
    @show(expr.head)#   :call
    @show(expr.args)#   Any[+, 1, 5]
    return (expr)    
end
ex = math_expr(+, 1, 5)
@show ((ex))#:((+)(1, 5))
@show eval(ex) #6 =#

##################################exemple2#################################
                #constructor not used in creating an expression: wrong use
#= function math_expr2(op, op1, op2)
    expr=:(:call, $op, $op1,$op2)  #expr lowercase is a var holds expression of list of things
    return (expr)    
end
ex2 = math_expr2(+, 1, 5)
@show ((ex2))#:((:call, +, 1, 5))
@show eval(ex2) #(:call, +, 1, 5)
@show eval(eval(ex2)) #(:call, +, 1, 5)
=#
##################################exemple22#################################
            #constructor not used in creating an expression: correct use
#= function math_expr22(op, op1, op2)
    expr=:(($op)($op1, $op2))  #expr lowercase is a var holds expression of list of things
    return (expr)    
end
ex22 = math_expr22(+, 1, 5)
@show ((ex22))# :((+)(1, 5))
@show eval(ex22) #6 =#

##################################exemple3#################################
 #constructor not used in creating an expression:  use of metaparse
#= function math_expr3(op, op1, op2)
   symbl=Symbol(op1,op,op2)
   # println((symbl))#1+5
   # println(typeof(symbl))#Symbol
    str=String(symbl)
   # display(str)# "1+5"
    expr=Meta.parse(str)
    display(expr);println()#:(1+5)
    return (expr)    
end
ex3 = math_expr3(+, 1, 5)
 @show (typeof(ex3))#:(1 + 5)
@show eval(ex2) #6 =#

####################################exemple4##################################""
#Base.:*(x::String, y::Int64) = string(x, y) # overload * for string*Int....string(String,Int) exists and concatenates  #@show "3"*2
#=
#nothing new but kept because it has an example of map and overloading * and string(String,int)
function math_expr2(op, op1, op2)
    (op1N,op2N)=map(x -> isa(x,Number) ? 2*x : x,(op1,op2))
    expr=Expr(:call, op, op1N, op2N)
    return expr    
end



ex2 = math_expr2(*, "e", 2)
@show (ex2)  # :((*)("e", 4))   which is the same as :("e"*4)
@show eval(ex2)# "e4"
 =#
