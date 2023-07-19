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
function extractJacDepNormal(varNum::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref # 
              push!(jacSet,  (a.args[2]))  # du[varNum=1]=rhs=u[5]+u[2] : 2 and 5 are stored in jacset
      end
      return a 
  end
  if length(jacSet)>0 jac[varNum]=jacSet end # jac={1->(2,5)}
end


function extractJacDepLoop(b::Int,niter::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  postwalk(rhs) do a   
      if a isa Expr && a.head == :ref # 
              push!(jacSet,  (a.args[2]))  #
      end
      return a
  end
  jac[:(($b,$niter))]=jacSet
end
