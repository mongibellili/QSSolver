
#= function modifyJacobian(jacobian :: SMatrix{N,N,Float64})where{N} #well test removing smatrix from begining 
  pos = zeros(SVector{N,SVector{N,Float64}})
  for i=1:N
    pos=setindex(pos,jacobian[i, :],i)
  end
  pos
end
function modifyJacobian2(jacobian :: SMatrix{N,N,Int})where{N} #well test removing smatrix from begining 
  pos = zeros(SVector{N,SVector{N,Int}})
  for i=1:N
    pos=setindex(pos,jacobian[i, :],i)  #@view
  end
  pos
end =#
#jacobian::SVector{T,SVector{T,Float64}} #  S-->D
# discreteJacobian::SVector{T,SVector{D,Float64}} #  d-->D
#ZC_jacobian::SVector{Z,SVector{T,Float64}}  #  S-->Z
#ZC_jacDiscrete::SVector{Z,SVector{D,Float64}} #  d-->Z
function createDependencyMatrix(jacobian::SVector{N,SVector{M,Float64}}) where{N,M}   # M effetcs N
   dep = zeros(SVector{M,SVector{N,Int}})
   for j=1:M
        arr=[]
        for i=1:N
          if jacobian[i][j]!=0#if jac[i,j] < -epselon ||  jac[i,j] > epselon # different than zero
                push!(arr,i)
          else
            push!(arr,0) # for some reason inner vectors want to have equal sizes... same when used tuples of tuples!!!
          end
        end
        dep=setindex(dep,arr,j)
    end  
    return dep  #dep=(tuple(arr...))  could use tuples; they also want equal size innner tuples!!!!!
end
function createDependencyToEventsDiscr(dD::SVector{D,SVector{T,Int}},dz::SVector{D,SVector{Z,Int}},eventHandlers::SVector{Y,EventHandlerStruct}) where{D,T,Z,Y}
  #we want Hd:  H (Y) effects d (numDiscr)
  #evDisc::SVector{D,Float64}
 # we want HZ1=Hd->dZ; [[..],[..],[..],[..]] svector{4,svector{2}}  Y=4  Z=2
  HZ1 = zeros(SVector{Y,SVector{Z,Int}})
  HD1 = zeros(SVector{Y,SVector{T,Int}})
 for j=1:Y
     evDiscr=eventHandlers[j].evDisc
    
     HZarr=[]
     HDarr=[]
     for i=1:D
            #------------event influences a discrete var
            if evDiscr[i]!=0# an event'j' effected a disc i
                  for k=1:Z
                    if length(HZarr)!=Z # first pass
                      if dz[i][k] !=0 #that disc 'i' effects ZC 'k'
                        push!(HZarr,k)
                      else
                        push!(HZarr,0) 
                      end
                    else # other passes...if you are here then the HZarr is already built with the wanted size
                      if dz[i][k] !=0 #that disc 'i' effects ZC 'k'
                        HZarr[k]=k # should i check if it is ==0 ie is it cheaper to reassign or if statment
                      end
                    end
                  end# end for k (ZC)
            end# end if event effects
    # end#end for i (discrete vars)
    # for i=1:D
           #------------event influences a cont var
           if evDiscr[i]!=0# an event'j' effected a disc i
            for k=1:T
              if length(HDarr)!=T # first pass
                if dD[i][k] !=0 #that disc 'i' effects der 'k'
                  push!(HDarr,k)
                else
                  push!(HDarr,0) 
                end
              else # other passes...if you are here then the HZarr is already built with the wanted size
                if dD[i][k] !=0 #that disc 'i' effects ZC 'k'
                  HDarr[k]=k # should i check if it is ==0 ie is it cheaper to reassign or if statment
                end
              end
            end# end for k (ZC)
           end# end if event effects
    end#end for i (cont vars)

     if length(HZarr)!=0
        HZ1=setindex(HZ1,HZarr,j)
     end
     if length(HDarr)!=0
         HD1=setindex(HD1,HDarr,j)
     end
 end #end for j  (events)
 return (HZ1,HD1)
   
end




function createDependencyToEventsCont(SD::SVector{T,SVector{T,Int}},sZ::SVector{T,SVector{Z,Int}},eventHandlers::SVector{Y,EventHandlerStruct}) where{T,Z,Y}
  #we want Hd:  H (Y) effects s (numContin)
  #evCont::SVector{D,Float64}
 # we want HZ1=Hd->sZ; [[..],[..],[..],[..]] svector{4,svector{2}}  Y=4  Z=2
  HZ2 = zeros(SVector{Y,SVector{Z,Int}})
  HD2 = zeros(SVector{Y,SVector{T,Int}})
 for j=1:Y
     evContin=eventHandlers[j].evCont
    
     HZarr=[]
     HDarr=[]
     for i=1:T
            #------------event influences a Continete var
            if evContin[i]!=0# an event'j' effected a Cont i
                  for k=1:Z
                    if length(HZarr)!=Z # first pass
                      if sZ[i][k] !=0 #that Cont 'i' effects ZC 'k'
                        push!(HZarr,k)
                      else
                        push!(HZarr,0) 
                      end
                    else # other passes...if you are here then the HZarr is already built with the wanted size
                      if sZ[i][k] !=0 #that Cont 'i' effects ZC 'k'
                        HZarr[k]=k # should i check if it is ==0 ie is it cheaper to reassign or if statment
                      end
                    end
                  end# end for k (ZC)
            end# end if event effects
    # end#end for i (Continete vars)
    # for i=1:D
           #------------event influences a cont var
           if evContin[i]!=0# an event'j' effected a Cont i
            for k=1:T
              if length(HDarr)!=T # first pass
                if SD[i][k] !=0 #that Cont 'i' effects der 'k'
                  push!(HDarr,k)
                else
                  push!(HDarr,0) 
                end
              else # other passes...if you are here then the HZarr is already built with the wanted size
                if SD[i][k] !=0 #that Cont 'i' effects ZC 'k'
                  HDarr[k]=k # should i check if it is ==0 ie is it cheaper to reassign or if statment
                end
              end
            end# end for k (ZC)
           end# end if event effects
    end#end for i (cont vars)

    if length(HZarr)!=0  #if no dependency leave initial svector intact
     HZ2=setindex(HZ2,HZarr,j)
    end
     if length(HDarr)!=0
     HD2=setindex(HD2,HDarr,j)
     end
 end #end for j  (events)
 return (HZ2,HD2)
   
end





function unionDependency(HZ1::SVector{Y,SVector{Z,Int}},HZ2::SVector{Y,SVector{Z,Int}})where{Z,Y}
  dep=zeros(SVector{Y,SVector{Z,Int}})
  for j=1:Y
          arr=[]
          for i=1:Z
                  if HZ1[j][i]!=0 ||  HZ2[j][i]!=0
                    push!(arr,i)
                  else
                    push!(arr,0)
                  end
          end
          dep=setindex(dep,arr,j)
  end  
  return dep

end




#= function creatediscreVarsDependencyToEvents(numDiscr::Int,eventHandlers::SVector{Y,EventHandlerStruct}) where{Y}
   #we want Hd:  H (Y) effects d (numDiscr)
   #evDisc::SVector{D,Float64}

   dep = zeros(SVector{Y,SVector{numDiscr,Int}})
  for j=1:Y
      evDiscr=eventHandlers[j].evDisc
      arr=[]
      for i=1:numDiscr
        if evDiscr[i]!=0#
              push!(arr,i)
        else
          push!(arr,0) # for some reason inner vectors want to have equal sizes... same when used tuples of tuples!!!
        end
      end
      dep=setindex(dep,arr,j)
  end  
  return dep
    
end

function createContVarsDependencyToEvents(numCont::Int,eventHandlers::SVector{Y,EventHandlerStruct}) where{Y}
  #we want H (Y) effects d (numCont)
  

  dep = zeros(SVector{Y,SVector{numCont,Int}})
 for j=1:Y
  evCont=eventHandlers[j].evCont
     arr=[]
     for i=1:numCont
       if evCont[i]!=0#
             push!(arr,i)
       else
         push!(arr,0) # for some reason inner vectors want to have equal sizes... same when used tuples of tuples!!!
       end
     end
     dep=setindex(dep,arr,j)
 end  
#=  for j=1:Y
  #evCont=eventHandlers[j].evCont
     #arr=[]
     for i=1:numCont
       if eventHandlers[j].evCont[i]!=0#
           eventHandlers[j].evCont=setindex(eventHandlers[j].evCont,i,i)
            #=  push!(arr,i)
       else
         push!(arr,0)  =#   #for some reason inner vectors want to have equal sizes... same when used tuples of tuples!!!
       end
     end
     dep=setindex(dep,eventHandlers[j].evCont,j)
 end  =#
 return dep 
   
end =#



function computeStates(initcond::SVector{N,Float64} ) where{N}
   return N#return length(initcond)   
end
function computeStates(initcond::MVector{N,Float64} ) where{N}
  return N  
end
function derivateJacobian(jacobian :: SMatrix{N,N,Float64})where{N}
return jacobian*jacobian #this is for  time-independent jacobians# run once: not needed cuz i used derQ instead
end

getOrderfromSolverMethod(::Val{1})=1
getOrderfromSolverMethod(::Val{2})=2
getOrderfromSolverMethod(::Val{3})=3
#getOrderfromSolverMethod(::Val{4})=1 #liqss1 #having more than 3 methods is causing an error 
#getOrderfromSolverMethod(::Val{5})=2  #liqss2