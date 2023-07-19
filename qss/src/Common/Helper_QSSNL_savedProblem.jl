



# these 2 function handle continuous problems only
function extractJacDepNormal(varNum::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},SD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref # 
              push!(jacSet,  (a.args[2]))  # du[varNum=1]=rhs=u[5]+u[2] : 2 and 5 are stored in jacset
              SDset=Set{Union{Int,Symbol,Expr}}()
              if haskey(SD, (a.args[2])) # if rhs var (5 for exple) already affects another var and it was visited before
                  SDset=get(SD,(a.args[2]),SDset) # get the set of vars that 5 affects
              end
              push!(SDset,  varNum) # 5 now affects 1 ...sdset=(...1)
              SD[(a.args[2])]=SDset  #SD={5->(...1)}
      end
      return a 
  end
  if length(jacSet)>0 jac[varNum]=jacSet end # jac={1->(2,5)}
end


function extractJacDepLoop(b::Int,niter::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},SD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  sdSet=Set{Union{Int,Symbol,Expr}}()
  postwalk(rhs) do a   
      if a isa Expr && a.head == :ref # 
              push!(jacSet,  (a.args[2]))  #
              if a.args[2] isa Int
                  SDset=Set{Union{Int,Symbol,Expr}}()
                  if haskey(SD, (a.args[2]))
                      SDset=get(SD,(a.args[2]),SDset)
                  end
                  push!(SDset,  :(($b,$niter)))
                  SD[(a.args[2])]=SDset
              elseif a.args[2] isa Symbol
                  push!(sdSet,  (a.args[2]))
              elseif a.args[2] isa Expr# The idea here will be explained in an example #for i=2:9 d[i]=u[i-1] end . 1 influences 2...In general I say if 2<=i+1<=9 ,i influences i+1: therefore + && - switch and * && / switch positions
                  temp=deepcopy(a.args[2])
                  if a.args[2].args[1]== :+
                      temp.args[1]= :-
                  elseif a.args[2].args[1]== :-
                      temp.args[1]= :+
                  elseif a.args[2].args[1]== :*
                      temp.args[1]= :/
                  elseif a.args[2].args[1]== :/
                      temp.args[1]= :*
                  end
                  push!(sdSet,  temp)
              end
      end
      return a
  end
  jac[:(($b,$niter))]=jacSet
  SD[:(($b,$niter))]=sdSet  # it should be the opposite but it is just a temporary storage and later checking if key is a tuple is easier
  #= for i=b:niter
      jac[i]=1
  end =#
end

# the following functions handle discrete problems
function extractJacDepNormal(varNum::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},jacDiscr :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},SD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  jacSet=Set{Union{Int,Symbol,Expr}}()
  jacDiscrSet=Set{Union{Int,Symbol,Expr}}()
  postwalk(rhs) do a   #
      if a isa Expr && a.head == :ref && a.args[1]==:q# 
          push!(jacSet,  (a.args[2]))  #
          SDset=Set{Union{Int,Symbol,Expr}}()
          if haskey(SD, (a.args[2]))
              SDset=get(SD,(a.args[2]),SDset)
          end
          push!(SDset,  varNum)
          SD[(a.args[2])]=SDset
      elseif a isa Expr && a.head == :ref && a.args[1]==:d# 
          push!(jacDiscrSet,  (a.args[2]))  #
          dDset=Set{Union{Int,Symbol,Expr}}()
          if haskey(dD, (a.args[2]))
              dDset=get(dD,(a.args[2]),dDset)
          end
          push!(dDset,  varNum)
          dD[(a.args[2])]=dDset
      end
      return a 
  end
  if length(jacSet)>0 jac[varNum]=jacSet end
  if length(jacDiscrSet)>0 jacDiscr[varNum]=jacDiscrSet end
end



function extractJacDepLoop(b::Int,niter::Int,rhs::Union{Int,Expr},jac :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},jacDiscr :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},SD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},dD :: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}) 
  #dD is for the saving-function case
  jacSet=Set{Union{Int,Symbol,Expr}}()
  jacDiscrSet=Set{Union{Int,Symbol,Expr}}()
  sdSet=Set{Union{Int,Symbol,Expr}}() #  when the index is a symbol or expr SD will be filled like Jac and later (the caller outside) will extract SDFunction. SDset uppercase is for numbers
  #dDSet=Set{Union{Int,Symbol,Expr}}()  # feature not implemented: I will restrict d[---] to integers ie  --- == intger
  postwalk(rhs) do a   
      if a isa Expr && a.head == :ref && a.args[1]==:q# 
              push!(jacSet,  (a.args[2]))  #
              if a.args[2] isa Int   #the index is a number
                  SDset=Set{Union{Int,Symbol,Expr}}()  #SDset   SD uppercase
                  if haskey(SD, (a.args[2]))
                      SDset=get(SD,(a.args[2]),SDset)
                  end
                  push!(SDset,  :(($b,$niter)))
                  SD[(a.args[2])]=SDset
              elseif a.args[2] isa Symbol   # there is only when symbol in this rhs
                  push!(sdSet,  (a.args[2]))
              elseif a.args[2] isa Expr   
                  temp=deepcopy(a.args[2])
                  if a.args[2].args[1]== :+
                      temp.args[1]= :-
                  elseif a.args[2].args[1]== :-
                      temp.args[1]= :+
                  elseif a.args[2].args[1]== :*
                      temp.args[1]= :/
                  elseif a.args[2].args[1]== :/
                      temp.args[1]= :*
                  end
                  push!(sdSet,  temp)
              end
      elseif a isa Expr && a.head == :ref && a.args[1]==:d
        push!(jacDiscrSet,  (a.args[2])) 
        if a.args[2] isa Int  #for now allow only d[integer]
          dDset=Set{Union{Int,Symbol,Expr}}()
          if haskey(dD, (a.args[2]))
              dDset=get(dD,(a.args[2]),dDset)
          end
          push!(dDset,  :(($b,$niter)))
          dD[(a.args[2])]=dDset
        #= elseif a.args[2] isa Symbol
            push!(dDSet,  (a.args[2]))
        elseif a.args[2] isa Expr
            temp=deepcopy(a.args[2])
            if a.args[2].args[1]== :+
                temp.args[1]= :-
            elseif a.args[2].args[1]== :-
                temp.args[1]= :+
            elseif a.args[2].args[1]== :*
                temp.args[1]= :/
            elseif a.args[2].args[1]== :/
                temp.args[1]= :*
            end
            push!(dDSet,  temp)
            =#
        end 
      end
      return a
  end
  jac[:(($b,$niter))]=jacSet
  jacDiscr[:(($b,$niter))]=jacDiscrSet
  SD[:(($b,$niter))]=sdSet  #symbol and expressions stored like jac

end










