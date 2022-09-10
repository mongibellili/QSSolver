
using MacroTools: postwalk
using TaylorSeries
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
  #input: expr=    :(x[$i-2]=$(odeExprs.args[i].args[2]))
  #output: x[$i-2][l=1]=$(odeExprs.args[i].args[2])[r=0]   # not using the coeffs methods i am using direct access with indices (shift -1)
  #       x[$i-2][l=2]=$(odeExprs.args[i].args[2])[r=1]/2  #second derivative
  
  # next task : it would be better if we can come up with a macro like @. bcast fusion to add coeffs[0] of all rhs then coeffs[1] ...
  #currently rhs=taylorA+taylorB+3+TaylorC*TaylorD is executed as a tree: rhs=taylorAB+3+taylorCD=taylorAB3+taylorCD=taylorAB3CD
  #which allocates in the middle because intermediary terms are taylor (heap-allocated) unlike floats (intermediary ops are stack-allocated)
  macro tayArth2(expr)# expr=    :(x[$i-2]=$(odeExprs.args[i].args[2]))
    Base.remove_linenums!(expr)
    lhs=postwalk(x -> x isa Expr && x.head == :ref && x.args[1] != :d ? getcoeffsLhs(x) : x, expr.args[1].args[1])#  !=d and head not needed 
   # @show lhs
    rhs=postwalk(x -> x isa Expr && x.head == :ref && x.args[1] != :d ? getcoeffsRhs(x) : x, expr.args[1].args[2])
    
    #@show rhs
    computedExpr=quote
                    local l=1  #local to avoid var leaking
                    local r=0
                    $lhs=$rhs
                    l=2
                    r=1
                    $lhs=$rhs/2
                           #  $lhsrhs
                            # println($lhsrhs)
                          #  println("runtime lhs= ",$lhs)
                end
    esc(computedExpr)
  end

  x = [Taylor1([0.0, 0.0,0.0]),Taylor1([0.0, 0.0,0.0])]
  q = [Taylor1([1.0, 1.0,1.5]),Taylor1([2.0, 2.0,2.5])]
  d = [1.0]
  display(@macroexpand @tayArth2  begin x[1]=9.8 -(d[1]/1.0)*(1e+5 *q[1]+30.0*q[2]) end)
   # @tayArth2  begin x[1]=(d[1])*(q[1]+3.0*q[2]) end
        @show x