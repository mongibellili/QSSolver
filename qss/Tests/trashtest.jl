using SymEngine
using MacroTools: postwalk
function eliminateRef(a)
  if a.args[2] isa Expr 
    if a.args[2].args[1]==:+
      a=Symbol((a.args[1]),(a.args[2].args[2]), "plus",(a.args[2].args[3]))
    elseif a.args[2].args[1]==:-
      a=Symbol((a.args[1]),(a.args[2].args[2]), "minus",(a.args[2].args[3]))
    end
  else
    a=Symbol((a.args[1]),(a.args[2]))
  end
  return a
end
function eliminateRef2(refEx)
  if refEx isa Expr 
    if refEx.args[1]==:+
      refEx=Symbol("q",(refEx.args[2]), "plus",(refEx.args[3]))
    elseif refEx.args[1]==:-
      refEx=Symbol("q",(refEx.args[2]), "minus",(refEx.args[3]))
    end
  else
    refEx=Symbol("q",(refEx))
  end
  return refEx
end

function changeVarNames_params(ex::Expr,symDict::Dict{Symbol,Expr})######maybe x is better for zc...if thats the case use this inside
  newEx=postwalk(ex) do element#postwalk to change var names and parameters
      if element isa Expr && element.head == :ref # 
        symarg=eliminateRef2(element.args[2])
        symDict[symarg]=element
      end
      return element
    end#end postwalk
  newEx
end

function restoreRef(coefExpr,symDict)
  newEx=postwalk(coefExpr) do element#postwalk to change var names and parameters
    if element isa Symbol && !(element in (:+,:-,:*,:/)) && haskey(symDict, element) 
      element=symDict[element]
     
    end
    return element
  end#end postwalk
  newEx
 
end

function test(symarg,rhsex,symDict::Dict{Symbol,Expr})
  m=postwalk(rhsex) do a   #after equs constructed, eliminate ref ; use new expr m since we still need z for equs below (below the caller of this)...postwalk
    if a isa Expr && a.head == :ref # a is expre of type A[n]  ie a=A[n]
      a=eliminateRef(a)
    end
    return a
  end
  #@show m
  symarg=eliminateRef2(symarg)
  #@show symarg
  basi = convert(Basic, m)
  #@show basi
  #createdSym=symbols("u$(a.args[2])")
  #@show createdSym
  coef = diff(basi, symarg)
  coefstr=string(coef)
  coefExpr=Meta.parse(coefstr)
  #dump(coefExpr)
  jacEntry=restoreRef(coefExpr,symDict)
  @show jacEntry
 @show rhsex
 #=  @show typeof(coef)
  @show (coef==112) =#
end
jac = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()
jac[1]=Set{Int}([1,2]);
symarg=:(i-1)
ex=:(2*q[i+1]*q[i-1]+5q[i-1]+q[i]+q[1])
symDict=Dict{Symbol,Expr}()
changeVarNames_params(ex,symDict)
#@show symDict
test(symarg,ex,symDict)

#= 
symarg=:(i-1)
ex=:(q[i-1])
test(symarg,ex)


symarg=:(1)
ex=:(q[1])
test(symarg,ex)


symarg=:(i)
ex=:(q[i])
test(symarg,ex)
 =#

