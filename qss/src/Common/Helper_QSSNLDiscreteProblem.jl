
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



function extractZCJacDepNormal(counter::Int,zcf::Expr,zcjac :: Vector{Vector{Int}},zcjacDiscr :: Vector{Vector{Int}},SZ ::Dict{Int,Set{Int}},dZ :: Dict{Int,Set{Int}}) 
  zcjacSet=Set{Int}()
  zcjacDiscrSet=Set{Int}()
  postwalk(zcf) do a   #
      if a isa Expr && a.head == :ref && a.args[1]==:q# 
          push!(zcjacSet,  (a.args[2]))  #
          SZset=Set{Int}()
          if haskey(SZ, (a.args[2]))
              SZset=get(SZ,(a.args[2]),SZset)
          end
          push!(SZset,  counter)
          SZ[(a.args[2])]=SZset
      elseif a isa Expr && a.head == :ref && a.args[1]==:d# 
          push!(zcjacDiscrSet,  (a.args[2]))  #
          dZset=Set{Int}()
          if haskey(dZ, (a.args[2]))
              dZset=get(dZ,(a.args[2]),dZset)
          end
          push!(dZset,  counter)
          dZ[(a.args[2])]=dZset
      end
      
      return a 
  end
  push!(zcjac,collect(zcjacSet))#convert set to vector
  push!(zcjacDiscr,collect(zcjacDiscrSet))
end








  function createDependencyToEventsDiscr(dD::Vector{Vector{Int}},dZ::Dict{Int64, Set{Int64}},eventDep::Vector{EventDependencyStruct}) 
    Y=length(eventDep)
lendD=length(dD)
    HD2 = Vector{Vector{Int}}(undef, Y)
    HZ2 = Vector{Vector{Int}}(undef, Y)
      for ii=1:Y
        HD2[ii] =Vector{Int}()# define it so i can push elements as i find them below
        HZ2[ii] =Vector{Int}()# define it so i can push elements as i find them below
      end
    for j=1:Y
      hdSet=Set{Int}()
      hzSet=Set{Int}()
        evdiscrete=eventDep[j].evDisc
        for i in evdiscrete
             if i<lendD  
              for k in dD[i]
                push!(hdSet,k)
              end
            end
              tempSet=Set{Int}()
              if haskey(dZ, i)
                tempSet=get(dZ,i,tempSet)
              end
              for kk in tempSet
                push!(hzSet,kk)
              end
        end
        HD2[j] =collect(hdSet)# define it so i can push elements as i find them below
        HZ2[j] =collect(hzSet)
    end #end for j  (events)
    return (HZ2,HD2)
  end 



function createDependencyToEventsCont(SD::Vector{Vector{Int}},sZ::Dict{Int64, Set{Int64}},eventDep::Vector{EventDependencyStruct}) 
  Y=length(eventDep)

  HD2 = Vector{Vector{Int}}(undef, Y)
  HZ2 = Vector{Vector{Int}}(undef, Y)
    for ii=1:Y
      HD2[ii] =Vector{Int}()# define it so i can push elements as i find them below
      HZ2[ii] =Vector{Int}()# define it so i can push elements as i find them below
    end
  for j=1:Y
    hdSet=Set{Int}()
    hzSet=Set{Int}()
      evContin=eventDep[j].evCont
      for i in evContin
            for k in SD[i]
              push!(hdSet,k)
            end
            tempSet=Set{Int}()
            if haskey(sZ, i)
              tempSet=get(sZ,i,tempSet)
            end
            for kk in tempSet
              push!(hzSet,kk)
            end
      end
      HD2[j] =collect(hdSet)# define it so i can push elements as i find them below
      HZ2[j] =collect(hzSet)
  end #end for j  (events)
  return (HZ2,HD2)
end 



#function unionDependency(HZ1::SVector{Y,SVector{Z,Int}},HZ2::SVector{Y,SVector{Z,Int}})where{Z,Y}
function unionDependency(HZ1::Vector{Vector{Int}},HZ2::Vector{Vector{Int}})
  Y=length(HZ1)


  HZ = Vector{Vector{Int}}(undef, Y)
    for ii=1:Y
      HZ[ii] =Vector{Int}()# define it so i can push elements as i find them below
    end
  for j=1:Y
    hzSet=Set{Int}()
    for kk in HZ1[j]
      push!(hzSet,kk)
    end
    for kk in HZ2[j]
      push!(hzSet,kk)
    end
    HZ[j]=collect(hzSet)
  end
  HZ 
end
