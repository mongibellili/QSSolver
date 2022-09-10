


 module module1
 using RuntimeGeneratedFunctions
 using StaticArrays
 using SymEngine
 using MacroTools: postwalk
 RuntimeGeneratedFunctions.init(@__MODULE__)

 
 struct EventHandlerStruct{T,D} #<: AbstractEventHandler
    id::Int
    evCont::SVector{T,Float64} 
    evDisc::SVector{D,Float64}
end



struct ODEProblem{T,D,Z,Y} 
    initConditions ::SVector{T,Float64}  
    discreteVars::MVector{D,Float64} 
   # jacobian::SVector{T,SVector{T,Float64}} 
    jacobian::SVector{T,SVector{T,Basic}} 
    eqs::Vector{Expr}
    discreteJacobian::SVector{T,SVector{D,Float64}} 
    #inputVars::SVector{T,Float64} 
    ZC_jacobian::SVector{Z,SVector{T,Float64}}  
    ZC_jacDiscrete::SVector{Z,SVector{D,Float64}} 
   # ZCinputVars::SVector{Z,Float64}   # in case input signal are function of t...use SVector{Z,Function} and create a different struct 
    eventHandlers::SVector{Y,EventHandlerStruct}
end
 
 



     macro ode(odeExprs)
         Base.remove_linenums!(odeExprs)   
         M=length(odeExprs.args)
         T=length((odeExprs.args[1].args[2]).args)  # number of cont vars.....correponds to the number of diff equ 
         D=length((odeExprs.args[2].args[2]).args)   # number of discrete vars

         exs = [quote $(Symbol(:f_, i))(u::Vector{Float64}) = $(odeExprs.args[i].args[2]) end for i = 3:2+T]

         odeExpr=postwalk(x -> x isa Expr && x.head==:ref  ? Symbol(:u,(x.args[2])) : x, odeExprs)
         println("-----modified-----=",odeExpr)

   #--------------------------------------------------------
                                    #               intial vals and diff equa part           #
                                    #---------------------------------------------------------
   
    contVars=SVector{T,Float64}((odeExpr.args[1].args[2]).args)
    discrVars=MVector{D,Float64}((odeExpr.args[2].args[2]).args)
    
    jac = zeros(SVector{T,SVector{T,Basic}}) 
    #jac = zeros(SVector{T,SVector{T,Float64}}) 
    jacDiscrete = zeros(SVector{T,SVector{D,Float64}}) 
   # inpuVarArr=[]
    u=[symbols("u$i") for i in 1:T] #number of symbols for cont vars
    d=[symbols("d$i") for i in 1:D] #number of symbols for cont vars
    println("vector of symbols= ",u)
    for i=3:2+T   # warning the number of diff equ might be less than T 
        jacArr=[]
        jacDiscArr=[]
        rhs=((odeExpr.args[i].args[2]))
        println("rhs= ",(rhs))
        basi = convert(Basic, rhs)#:(2u1+3u2))
        #extract jaco components
        for j=1:T            # the number of diffequ always coincides with the number of continuous vars?
            coef=diff(basi,u[j])
            #println(typeof(coef))
            push!(jacArr,coef)
        end
       #display(jacArr)
        #extract discreteJaco components
        for j=1:D            # problem here since we do not know the number of discrete vars
            coef=diff(basi,d[j])
            push!(jacDiscArr,coef)
        end
        #extract inpu vars
        #= for j=1:T
            basi=subs(basi,u[j]=>0)          
        end
        for j=1:D  # this could ve been done in the previous 1:D loop
            basi=subs(basi,d[j]=>0)
        end =#
        jac=setindex(jac,jacArr,i-2)
        jacDiscrete=setindex(jacDiscrete,jacDiscArr,i-2)
       # push!(inpuVarArr,basi) 
    end  
    #inputVars=SVector{T,Float64}(inpuVarArr)
   





    Z=M-2-T                         #number of if statement ===number of ZC
    ZC_jac = zeros(SVector{Z,SVector{T,Float64}}) 
    ZC_jacDiscrete = zeros(SVector{Z,SVector{D,Float64}}) 
    
    #ZC_inpuVarArr=[]            # temp helper arr
    evsArr=[]            # temp helper arr to gather event stucts
    u=[symbols("u$i") for i in 1:T] #number of symbols for cont vars
    d=[symbols("d$i") for i in 1:D] #number of symbols for disc vars
   
    for i=3+T:M   # loop that parses all if statments

                                    #----------------------------------------------------
                                    #                       Zero crossing part           #
                                    #-----------------------------------------------------
        jacArr=[]
        jacDiscArr=[]
        zc=odeExpr.args[i].args[1].args[2] 
        #dump(zc)
        basi = convert(Basic, zc)#:(2u1+3u2))
        #----------------------extract jaco components
        for j=1:T            # the number of odeProblem always coincides with the number of continuous vars
            coef=diff(basi,u[j])#SymEngine.diff(basi,u[j])
            push!(jacArr,coef)
        end
        ZC_jac=setindex(ZC_jac,jacArr,i-2-T)
        #------------------------extract discreteJaco components
        for j=1:D            # 
            coef=diff(basi,d[j])
            push!(jacDiscArr,coef)
        end
        ZC_jacDiscrete=setindex(ZC_jacDiscrete,jacDiscArr,i-2-T)
       #=  #--------------------------extract inpu vars
        for j=1:T
            basi=subs(basi,u[j]=>0)          
        end
        for j=1:D  # this could ve been done in the previous 1:D loop
            basi=subs(basi,d[j]=>0)
        end
        push!(ZC_inpuVarArr,basi) =#
                                    #----------------------------------------------------
                                    #                       events part                 #
                                    #-----------------------------------------------------
        #println("if statement number ",i-2)
        posEvExp=odeExpr.args[i].args[2]  # i corresponds to zc number i; it has 2 events (arg[2]=posEv and arg[3]=NegEv)
        negEvExp=odeExpr.args[i].args[3]
        k=i-2-T
        indexPosEv=2*k-1
        indexNegEv=2*k
        # now each ev can have many statements 
        #println(length(posEvExp.args))
                         #------------------pos Event--------------------#
        posEv_disArr=@SVector zeros(D) #  discrete
        posEv_conArr=@SVector zeros(T)  #  continous
        for j=1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
           posEvLHS=posEvExp.args[j].args[1]
           posEvRHS=posEvExp.args[j].args[2]
           #dump(posEvLHS)# symbol at the LHS
           basicLHS = convert(Basic, posEvLHS)
           
           discVarpositionArray = indexin(basicLHS, d)
           
           if !(discVarpositionArray[1] === nothing)
            posEv_disArr=setindex(posEv_disArr,posEvRHS,discVarpositionArray[1])
           else # lhs is not a disc var 
                conVarpositionArray = indexin(basicLHS, u)
                if !( conVarpositionArray[1]=== nothing)
                    posEv_conArr=setindex(posEv_conArr,posEvRHS,conVarpositionArray[1])
                else
                        println("LHS is neither a cont nor a discr var!!")
                end
            end
        end

              #------------------neg Event--------------------#
              negEv_disArr=@SVector zeros(D) #  discrete
              negEv_conArr=@SVector zeros(T)  #  continous
              for j=1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
                 negEvLHS=negEvExp.args[j].args[1]
                 negEvRHS=negEvExp.args[j].args[2]
                 #dump(negEvLHS)# symbol at the LHS
                 basicLHS = convert(Basic, negEvLHS)
                 
                 discVarpositionArray = indexin(basicLHS, d)
                 
                 if !(discVarpositionArray[1] === nothing)
                  negEv_disArr=setindex(negEv_disArr,negEvRHS,discVarpositionArray[1])
                 else # lhs is not a disc var 
                      conVarpositionArray = indexin(basicLHS, u)
                      if !( conVarpositionArray[1]=== nothing)
                          negEv_conArr=setindex(negEv_conArr,negEvRHS,conVarpositionArray[1])
                      else
                              println("LHS is neither a cont nor a discr var!!")
                      end
                  end
              end


       structposEvent=EventHandlerStruct(indexPosEv,posEv_conArr,posEv_disArr)
       push!(evsArr,structposEvent)
       structnegEvent=EventHandlerStruct(indexNegEv,negEv_conArr,negEv_disArr)
       push!(evsArr,structnegEvent)
     

    end  # end for that parses all if-statments 

     Y=2*Z
     eventHandlers=SVector{Y,EventHandlerStruct}(evsArr)  # 2*Z each zc yields 2 events
    #println(eventHandlers)
   # ZCinputVars=SVector{Z,Float64}(ZC_inpuVarArr)
    #println("------------------------ ",jac)
    #based on the type of the problem after a different user input call the appropriate struct
    myodeProblem=ODEProblem(contVars,discrVars,jac,exs,jacDiscrete,ZC_jac,ZC_jacDiscrete, eventHandlers)
    #myodeProblem=ODEProblem(contVars,discrVars,jac,exs,jacDiscrete,inputVars,ZC_jac,ZC_jacDiscrete,ZCinputVars, eventHandlers)

    myodeProblem





         
     end
     
 #=     function prepareProb()
        # display(eqs[1].args[2])
       
       arr
     end =#
     function integrate(odep::ODEProblem)


#*********************************problem info*****************************************
initConditions = odep.initConditions  
discVars=odep.discreteVars
#----------to compute derivatives
jacobian = odep.jacobian
discJac=odep.discreteJacobian
#inputVars=odep.inputVars
#----------to compute ZC expressions
zc_jac=odep.ZC_jacobian
ZC_discJac=odep.ZC_jacDiscrete
#ZC_input=odep.ZCinputVars
#-----------to execute event handlers
evHandlr=odep.eventHandlers
println("jacobian= ",jacobian)
println("discJac= ",discJac)
#println("inputVars= ",inputVars[2])
#=   t=Taylor1(2)
println("inputVars= ",inputVars[2](t)) =#
println("zc_jac= ",zc_jac)
println("ZC_discJac= ",ZC_discJac)
#println("ZC_input= ",ZC_input)
println("evHandlr= ",evHandlr)


numCont = computeStates(initConditions)
numDiscr = computeStates(discVars)
#numberZC = computeStates(ZC_input)
#=   println("numCont= ",numCont)
println("numDiscr= ",numDiscr)
println("numberZC= ",numberZC) =#
SD =createDependencyMatrix(jacobian)
#println("SD= ",SD)
#= dD =createDependencyMatrix(discJac) # temp dependency to be used to determine HD1 and HZ1 HD=Hd-dD Union Hs-sD
SZ =createDependencyMatrix(zc_jac) 
#println("SZ= ",SZ)
dZ =createDependencyMatrix(ZC_discJac) # temp dependency to be used to determine HD2 and HZ2
#=   println(SD)
println(dD)
println(sZ)
println(dZ) =#
#Hd=creatediscreVarsDependencyToEvents(numDiscr,evHandlr) #temp
#println(Hd)
#Hs=createContVarsDependencyToEvents(numCont,evHandlr)
#println(Hs)
HZ1HD1=createDependencyToEventsDiscr(dD,dZ,evHandlr) 
#=   println(HZ1HD1[1])
println(HZ1HD1[2]) =#
HZ2HD2=createDependencyToEventsCont(SD,SZ,evHandlr) 
#=   println(HZ2HD2[1])
println(HZ2HD2[2]) =#
#HD=Hd->dD Union Hs->sD
#HZ=Hd->dZ Union Hs->sZ

HZ=unionDependency(HZ1HD1[1],HZ2HD2[1])
#println("HZ= ",HZ)
HD=unionDependency(HZ1HD1[2],HZ2HD2[2])
 =#

println("SD= ",SD)
#= println("SZ= ",SZ)
println("HZ= ",HZ)
println("HD= ",HD) =#








        arr=Vector{Function}()
        for i=1:length(odep.eqs)
        # push!(arr,eval(eqs[i])) #f_1
         push!(arr, @RuntimeGeneratedFunction(odep.eqs[i].args[2]))
        end
        u=[2.2,3.4]
         for k=1:1
             for i=1:length(u)
                 derivate(i,u,arr[i])
             end
         end
         println(u)
        
     end
   function derivate(i::Int,u::Vector{Float64},f::Function)
         u[i]=f(u)
         return nothing
     end



     function modifyJacobian(jacobian :: SMatrix{N,N,Float64})where{N} #well test removing smatrix from begining 
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
      end
      #jacobian::SVector{T,SVector{T,Float64}} #  S-->D
      # discreteJacobian::SVector{T,SVector{D,Float64}} #  d-->D
      #ZC_jacobian::SVector{Z,SVector{T,Float64}}  #  S-->Z
      #ZC_jacDiscrete::SVector{Z,SVector{D,Float64}} #  d-->Z
      function createDependencyMatrix(jacobian::SVector{N,SVector{M,Basic}}) where{N,M}   # M effetcs N
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
      
      
      
      
      
      function computeStates(initcond::SVector{N,Float64} ) where{N}
         return N#return length(initcond)   
      end
      function computeStates(initcond::MVector{N,Float64} ) where{N}
        return N  
      end
      function derivateJacobian(jacobian :: SMatrix{N,N,Float64})where{N}
      return jacobian*jacobian #this is for  time-independent jacobians# run once: not needed cuz i used derQ instead
      end
      


 end#end module
 #-------------------------user space-----------------------------------
 using BenchmarkTools
 function test()
 odeprob=module1.@ode begin
    u=[2.2,3.4,2.2]
    d=[1.0,0.5]
     du[1]=u[1]+2.0u[2]*u[1]   
     du[2]=u[2]-3
     du[3]=u[3]+u[2]-u[1]+cos(t)
  #=    du1=u1+2.0u2   
     du2=u2-3#*cos(t) =#



#=     
     
     du1=u3-u2*u1+2.0u2     # du1...2...N are expected to be in order....later can fix this
     du2=-u1-3u2-d1
     du3=u3+u2+d1+d2   
     if u1+0.7+d2 >0   #still have not added 'usercode checking'....throw error msg
         d1=0.1
     else
         d1=1.0
         u2=4.5
     end
     if u2+d1 >0 
         d2=0.33
         u3=2.2
     else
         d2=1.0
     end =#

 end  
 
 #eqs=module1.prepareProb(odeprob)
 #@btime module1.integrate(u,eqs)
 module1.integrate(odeprob)
 end
 #@btime test()
 test()
 
 


