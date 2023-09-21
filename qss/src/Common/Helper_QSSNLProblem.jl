function changeExprToFirstValue(ex::Expr)##
  newEx=postwalk(ex) do a  # change u[1] to u[1][0]
      if a isa Expr && a.head == :ref && a.args[1]==:q
           outerRef=Expr(:ref)
          push!(outerRef.args,a)
          push!(outerRef.args,:(0))
          a=outerRef
      end
      return a
  end

  newEx
  end

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

function restoreRef(coefExpr,symDict)
  newEx=postwalk(coefExpr) do element#postwalk to change var names and parameters
    if element isa Symbol && !(element in (:+,:-,:*,:/)) && haskey(symDict, element) 
      element=symDict[element]
      element=changeExprToFirstValue(element)
    end
    return element
  end#end postwalk
  newEx
 
end
function changeVarNames_params(ex::Expr,stateVarName::Symbol,muteVar::Symbol,param::Dict{Symbol,Number},symDict::Dict{Symbol,Expr})######maybe x is better for zc...if thats the case use this inside
  newEx=postwalk(ex) do element#postwalk to change var names and parameters
      if element isa Symbol   
          if haskey(param, element)#symbol is a parameter
              element=param[element] 
          elseif element==stateVarName #symbol is a var
              element=:q 
          elseif element==:discrete #symbol is a discr var
              element=:d
          elseif element==muteVar #symbol is a mute var
              element=:i
          #= else  # + - * /
               =#
          end
      elseif element isa Expr && element.head == :ref # 
            symarg=eliminateRef2(element.args[2])
            symDict[symarg]=element
          
      end
      return element
    end#end postwalk
  newEx
end

function changeVarNames_params(ex::Expr,stateVarName::Symbol,muteVar::Symbol,param::Dict{Symbol,Number})######maybe x is better for zc...if thats the case use this inside
  newEx=postwalk(ex) do element#postwalk to change var names and parameters
      if element isa Symbol   
          if haskey(param, element)#symbol is a parameter
              element=param[element] 
          elseif element==stateVarName #symbol is a var
              element=:q 
          elseif element==:discrete #symbol is a discr var
              element=:d
          elseif element==muteVar #symbol is a mute var
              element=:i
          #= else  # + - * /
               =#
          end
      end
      return element
    end#end postwalk
  newEx
end

# these 2 function handle continuous problems only
function extractJacDepNormal(varNum::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}, exacteJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  m=postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref # 
              push!(jacSet,  (a.args[2]))  # du[varNum=1]=rhs=u[5]+u[2] : 2 and 5 are stored in jacset
              a=eliminateRef(a)
      end
      return a 
  end
  basi = convert(Basic, m)
  for i in jacSet
    symarg=eliminateRef2(i)
    
    #@show basi
    #createdSym=symbols("u$(a.args[2])")
    #@show createdSym
  #=   coef1 = diff(basi, symarg)#df
    coef2 = diff(coef1, symarg)#ddf
    coefstr=string(coef1,-,"(",coef2,")*",symarg,*,0.5) =#

    coef = diff(basi, symarg)
    coefstr=string(coef)

   # coefstr=string("(",basi,")/",symarg) 


   # @show coefstr
    
    coefExpr=Meta.parse(coefstr)
    #dump(coefExpr)
    jacEntry=restoreRef(coefExpr,symDict)
   # @show jacEntry
    exacteJacExpr[:(($varNum,$i))]=jacEntry
  end
  if length(jacSet)>0 jac[varNum]=jacSet end # jac={1->(2,5)}
end




function extractJacDepLoop(b::Int,niter::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}} ,exacteJacExpr :: Dict{Expr,Union{Float64,Int,Symbol,Expr}},symDict::Dict{Symbol,Expr}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  m=postwalk(rhs) do a   
      if a isa Expr && a.head == :ref # 
              push!(jacSet,  (a.args[2]))  #
              a=eliminateRef(a)
      end
      return a
  end
  basi = convert(Basic, m)
  for i in jacSet
    symarg=eliminateRef2(i);

    coef = diff(basi, symarg)
    coefstr=string(coef);

    #= coef1 = diff(basi, symarg)#df
    coef2 = diff(coef1, symarg)#ddf
    coefstr=string(coef1,-,"(",coef2,")*",symarg,*,0.5) =#

   
   # coefstr=string("(",basi,")/",symarg) 

    coefExpr=Meta.parse(coefstr)
    jacEntry=restoreRef(coefExpr,symDict)
    exacteJacExpr[:((($b,$niter),$i))]=jacEntry
  end
  if length(jacSet)>0 jac[:(($b,$niter))]=jacSet end
end