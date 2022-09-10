


 module module1
 using RuntimeGeneratedFunctions
 using StaticArrays
 using SymEngine
 using TaylorSeries

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

         exs = [quote $(Symbol(:f_, i))(q::Vector{Taylor1{Float64}},t::Taylor1{Float64}) = $(odeExprs.args[i].args[2]) end for i = 3:2+T]

         odeExpr=postwalk(x -> x isa Expr && x.head==:ref  ? Symbol(:q,(x.args[2])) : x, odeExprs)
        # println("-----modified-----=",odeExpr)

   #--------------------------------------------------------
                                    #               intial vals and diff equa part           #
                                    #---------------------------------------------------------
   
    contVars=SVector{T,Float64}((odeExpr.args[1].args[2]).args)
    discrVars=MVector{D,Float64}((odeExpr.args[2].args[2]).args)
    
    jac = zeros(SVector{T,SVector{T,Basic}}) 
    #jac = zeros(SVector{T,SVector{T,Float64}}) 
    jacDiscrete = zeros(SVector{T,SVector{D,Float64}}) 
   # inpuVarArr=[]
    u=[symbols("q$i") for i in 1:T] #number of symbols for cont vars
    d=[symbols("d$i") for i in 1:D] #number of symbols for cont vars
    #println("vector of symbols= ",u)
    for i=3:2+T   # warning the number of diff equ might be less than T 
        jacArr=[]
        jacDiscArr=[]
        rhs=((odeExpr.args[i].args[2]))
       # println("rhs= ",(rhs))
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
    u=[symbols("q$i") for i in 1:T] #number of symbols for cont vars
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
            coef=diff(basi,u[j])
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
#= println("jacobian= ",jacobian)
println("discJac= ",discJac)
#println("inputVars= ",inputVars[2])
#=   t=Taylor1(2)
println("inputVars= ",inputVars[2](t)) =#
println("zc_jac= ",zc_jac)
println("ZC_discJac= ",ZC_discJac)
#println("ZC_input= ",ZC_input)
println("evHandlr= ",evHandlr) =#


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

#println("SD= ",SD)
#= println("SZ= ",SZ)
println("HZ= ",HZ)
println("HD= ",HD) =#








        f=Vector{Function}()
        for i=1:length(odep.eqs)
        # push!(arr,eval(eqs[i])) #f_1
         push!(f, @RuntimeGeneratedFunction(odep.eqs[i].args[2]))
        end
       # initCond=[2.2,3.4]
         
        
     
     order=2
     Debug=false
     # Debug=true
       ft=10.0
       absQ=1e-6
       relQ=1e-3
       dep= @SVector[[0,2],[1,2]]
       #dep= @SVector[[0]]
       states=length(initConditions)
       x=Taylor1{Float64}[]
       q=Taylor1{Float64}[]
       t=Taylor1(order)
       #t=Taylor1([1e-120,1,0],order-1) # this is bad cuz 1/t explodes
     #  initTime=5.562684646268003e-308 ########**********sqrt and xreal exponents complain about expansion around zero....this is bad cuz solution explodes (tested)
     #  initTime= 1e-152 # only a not too small number works
       initTime=0.0
       interpolationOrder=1 # 1 means savedVars = x + t derX
       arr = []
       for i = 1:states
            push!(arr, [])
       end
       savedVars = SVector{states,Array{Taylor1{Float64}}}(tuple(arr...))    
       savedTimes = Array{Float64}([initTime])
       tx = @MVector zeros(states)
       tq = @MVector zeros(states)
       nextStateTime = @MVector zeros(states)
       nextInputTime = @MVector zeros(states)
       quantum = @MVector zeros(states)
       for i=1:states
           push!(x,Taylor1(zeros(order+1),order)+initConditions[i])
           push!(q,Taylor1(zeros(order),order-1))
           push!(savedVars[i],Taylor1(zeros(interpolationOrder+1),interpolationOrder))
       end
       if Debug println("initialsavedVars  = ",savedVars) end
        for k=1:order # compute initial derivatives for x and q
          # println("k= ",k)
           for i=1:states  
               q[i].coeffs[k] =x[i].coeffs[k]
           end
         # println("k $k q= ",q) 
          for i=1:states 
           x[i].coeffs[k+1] =((differentiate(f[i](q,t+initTime),k-1)).coeffs[1] )/factorial(k) # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
          end
       end     
       if Debug       
       println("initialX with derivatives= ",x)
       println("initialq with derivatives= ",q)
       end
      
       for i = 1:states
           tx[i] = initTime
           tq[i] = initTime
           for j=1:interpolationOrder
               savedVars[i][1].coeffs[j]=x[i].coeffs[j]
           end
           quantum[i] = relQ * abs(x[i].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)
           
           if quantum[i] < absQ
               quantum[i] = absQ
           end
           
           computeNextTime(Val(2), i, initTime, nextStateTime, x, quantum)
           computeNextInputTime(Val(2), i, initTime,0.0,0.0+t,f[i], nextInputTime, x,q, quantum)
       end

       if Debug 
         println("initial quantum= ",quantum)
         println("intial nextstateTime= ",nextStateTime) 
         println("intial nextInputTime= ",nextInputTime) 
       end
       simt = initTime
       count=0
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
       while simt < ft && count < 15555
          count+=1

          sch = updateScheduler(nextStateTime,nextInputTime)
          simt= sch[2]
           index = sch[1]

           #  simt, index = findmin(nextStateTime)
           if Debug 
             println("loop$count : nextStateTime rightafter findmin = ",nextStateTime)
             println("loop$count : nextInputTime rightafter findmin = ",nextInputTime)
             println("loop$count :  ",sch[3])
            end
           if Debug  println("simtime= ",simt," index=  ",index) end           
     # println(x[index](t))
             elapsed = simt - tx[index]
   #  if sch[3] ==:ST_STATE

          
            #shift taylor
           for k=1:order #evaluate x at new time, and derivatives...remember coeffs store the over facto also              
                x[index].coeffs[k] =(((differentiate(x[index],k-1)(t+elapsed))).coeffs[1])/factorial(k-1)
           end
          # x[index].coeffs[1]=x[index](t+elapsed).coeffs[1]
          # x[index].coeffs[2]=((differentiate(x[index]))(t+elapsed)).coeffs[1]
         
          
           if Debug println("loop $count :x after evaluatex = ",x) end
           quantum[index] = relQ * abs(x[index].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
           if quantum[index] < absQ
               quantum[index] = absQ
           end
           if Debug  println("loop $count :quantum = ",quantum) end
           tx[index]=simt
           for k=1:order 
               q[index].coeffs[k] =x[index].coeffs[k] #updateQ
           end
           if Debug println("loop $count :q after updateQ = ",q) end
           tq[index]=simt #tq needed for higher orders
           if sch[3] ==:ST_STATE
                       computeNextTime(Val(2), index, simt, nextStateTime, x, quantum) #
                         # computeNextInputTime(Val(2), index, simt,elapsed,elapsed+t,f[index], nextInputTime, x,q, quantum)
                         if Debug println("loop $count : nextStateTime after computenexttime= ",nextStateTime) end
                               for i = 1:length(dep[index])
                             
                               j = dep[index][i] # this line before the next line or vice versa gave the same bench results
                               if j != 0   
                                 # println("dependency loop $j= ")              
                                 elapsed = simt - tx[j]                 
                                 if elapsed > 0
                                   #@timeit "integrate state" 
                                   x[j].coeffs[1]=x[j](t+elapsed).coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
                                   q[j].coeffs[1] =x[j].coeffs[1]
                                   tx[j] = simt

                                 end
                                 if Debug 
                                   println("loop $count : x before derivative = ",x) 
                                 println("loop $count : q before derivative = ",q) 
                                 end
                                 #@timeit "compute deriv" 
                                 computeDerivative(Val(2),j, f[j], x, q, simt+t)
                                 #@timeit "Recompute next" 
                                 if Debug 
                                   println("loop $count : x after derivative = ",x) 
                                   println("loop $count : q after derivative = ",q) 
                                 end
                                 reComputeNextTime(Val(2),j, simt, nextStateTime, x, q, quantum)
                                 if Debug println("loop $count : nextStateTime after reComputenexttime= ",nextStateTime) end
                             
                               end#end if j!=0
                         end#end for
           else  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself
                       # computeNextTime(Val(2), index, simt, nextStateTime, x, quantum) #
                      # computeDerivative(Val(2),index, f[index], x, q, simt+t)
                       computeNextInputTime(Val(2), index, simt,elapsed,simt+t,f[index], nextInputTime, x,q, quantum)
                         if Debug println("under inputchange loop $count : nextInputTime after computenextInputtime= ",nextInputTime) end
                         computeDerivative(Val(2),index, f[index], x, q, simt+t)
                           
                           #@timeit "Recompute next" 
                           if Debug 
                             println("under inputchange loop $count : x after derivative = ",x) 
                             println("under input change loop $count : q after derivative = ",q) # q not updated here
                           end
                          # reComputeNextTime(Val(2),index, simt, nextStateTime, x, q, quantum)
                           for i = 1:length(dep[index])
                             
                               j = dep[index][i] # this line before the next line or vice versa gave the same bench results
                               if j != 0   
                               # println("dependency loop $j= ")              
                                 elapsed = simt - tx[j]                 
                                 if elapsed > 0
                                   #@timeit "integrate state" 
                                   x[j].coeffs[1]=x[j](t+elapsed).coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
                                   q[j].coeffs[1] =x[j].coeffs[1]
                                   tx[j] = simt
               
                                 end
                                 if Debug 
                                   println("under inputchange loop $count : x before derivative = ",x) 
                                 println("under inputchange loop $count : q before derivative = ",q) 
                                 end
                                 #@timeit "compute deriv" 
                                 computeDerivative(Val(2),j, f[j], x, q, simt+t)
                                 #@timeit "Recompute next" 
                                 if Debug 
                                   println("under inputchange loop $count : x after derivative = ",x) 
                                   println("under inputchange loop $count : q after derivative = ",q) 
                                 end
                                 reComputeNextTime(Val(2),j, simt, nextStateTime, x, q, quantum)
                                 if Debug println("under inputchange loop $count : nextStateTime after reComputenexttime= ",nextStateTime) end
                             
                               end#end if j!=0
                           end#end for
               
           end

#=           else#input 
                 computeNextInputTime(Val(2), index, simt,elapsed,elapsed+t,f[index], nextInputTime, x,q, quantum)
                #=  for k=1:order #evaluate x at new time, and derivatives...remember coeffs store the over facto also              
                   x[index].coeffs[k] =(((differentiate(x[index],k-1)(t+elapsed))).coeffs[1])/factorial(k-1)
                 end =#
                 computeDerivative(Val(2),index, f[index], x, q, simt+t)
         end#end if state/input =#
           for k = 1:states
               temp=Taylor1(interpolationOrder)
               for j=1:interpolationOrder
                   temp.coeffs[j]=x[k].coeffs[j]
               end
              # println("temp= ",temp)
               push!(savedVars[k], temp)
           end
           push!(savedTimes, simt)
       end#end while
       (savedTimes, savedVars)
   end#end integrate





 



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
      
      function computeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}}, quantum::MVector{T,Float64})where{T}
        absDeltaT=1e-12
          if (x[i].coeffs[3]) != 0
              #= println("currentTime = ",currentTime)
              println("quantum[$i] = ",quantum[i])
              println("(x[$i].coeffs[3])*2 = ",(x[i].coeffs[3])*2) =#
              #nextTime[i] = currentTime + sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
              tempTime=max(sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))),absDeltaT)
              if tempTime!=absDeltaT #normal
                  nextTime[i] = currentTime + tempTime#sqrt(abs(quantum[i] / ((x[i].coeffs[3])*2))) #*2 cuz coeff contains fact()
              else#usual sqrt(quant/der) is very small
                #println("(x[$i].coeffs[3])*2 = ",(x[i].coeffs[3])*2)
                x[i].coeffs[3]=sign(x[i].coeffs[3])*(abs(quantum[i])/(absDeltaT*absDeltaT))/2# adjust second derivative if it is too high
                nextTime[i] = currentTime + tempTime
               # println("corrected (x[$i].coeffs[3])*2 = ",(x[i].coeffs[3])*2)
              end
             #println("schedule state at= ",nextTime[i])
          else
            nextTime[i] = Inf
          end
      end
      function computeNextInputTime(::Val{2}, i::Int, currentTime::Float64,elapsed::Float64, t::Taylor1{Float64} ,f::Function,nextInputTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}}, quantum::MVector{T,Float64})where{T}
        df=0.0
        oldDerDerX=((x[i].coeffs[3])*2)
        newDerDerX=((differentiate(f(q,t),1)).coeffs[1] )# do not put /factorial(2) cuz here we actually need derder not store the coeff
       # println("-------------------------------------------",f(q,t))
          if elapsed > 0
           # println("elapsed= ",elapsed) 
           # println("oldDerDerX= ",oldDerDerX) 
           # println("newDerDerX= ",newDerDerX) 
            df=(newDerDerX-oldDerDerX)/(elapsed/2)
            #df=newDerDerX
           # println("df= ",df) 
          else
            df= quantum[i]*1e12
          end
            
         if df!=0
          nextInputTime[i]=currentTime+(abs(quantum[i] / df))^(1/3)
          
         else
          nextInputTime[i] = Inf
          end
  #=         oldDerX=((x[i].coeffs[2]))
          newDerX=((f(q,t)).coeffs[1] )
            if elapsed > 0
              df=(newDerX-oldDerX)/elapsed
              println("df= ",df) 
            else
              df= quantum[i]*1e6
            end
              
           if df!=0
            nextInputTime[i]=currentTime+(abs(quantum[i] / df))^(1/2)
           else
            nextInputTime[i] = Inf
            end =#
      end
      function computeDerivative( ::Val{2} ,j::Int, f::Function ,x::Vector{Taylor1{Float64}} , q::Vector{Taylor1{Float64}}, t::Taylor1{Float64}  )
          #x[j].coeffs[2] =((differentiate(f(q,t),0)).coeffs[1] )/factorial(1) 
          x[j].coeffs[2] =(f(q,t)).coeffs[1]
          x[j].coeffs[3] =((differentiate(f(q,t),1)).coeffs[1] )/factorial(2)
          #println("x[j].coeffs[2]= ",x[j].coeffs[2])
          return nothing
      end
      function reComputeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::Vector{Taylor1{Float64}},q::Vector{Taylor1{Float64}}, quantum::MVector{T,Float64})where{T}
          coef=@SVector [q[index].coeffs[1] - (x[index].coeffs[1]) - quantum[index], q[index].coeffs[2]-x[index].coeffs[2],-(x[index].coeffs[3])*2]
         #=  println("coef[1]= ",coef[1])
          println("coef[2]= ",coef[2])
          println("coef[3]= ",coef[3]) =#
          time1 = currentTime + minPosRoot(coef, Val(2))
         # println("time1= ",time1)
          coef=setindex(coef,q[index].coeffs[1] - (x[index].coeffs[1]) + quantum[index],1)
         # println("coef[1]= ",coef[1])
          time2 = currentTime + minPosRoot(coef, Val(2))
        #  println("time2= ",time2)
          nextTime[index] = time1 < time2 ? time1 : time2
          return nothing
      end
      function minPosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
          mpr=-1 #coef1=c, coef2=b, coef3=a
          if coeff[3] == 0 #|| (10000 * abs(coeff[3])) < abs(coeff[2])
              if coeff[2] == 0
                mpr = Inf
              else 
                mpr = -coeff[1] / coeff[2]
              end
              if mpr < 0
                mpr = Inf
              end
          else 
             #double disc;
              disc = coeff[2] * coeff[2] - 4 * coeff[3] * coeff[1]#b^2-4ac
              if disc < 0 # no real roots
                mpr = Inf
              else 
                #double sd, r1;
                sd = sqrt(disc);
                r1 = (-coeff[2] + sd) / (2 * coeff[3]);
                if r1 > 0 
                  mpr = r1;
                else 
                  mpr = Inf;
                end
                r1 = (-coeff[2] - sd) / (2 * coeff[3]);
                if ((r1 > 0) && (r1 < mpr)) 
                  mpr = r1;
                end
              end
              
          end
          return mpr
      end
  
      function updateScheduler(nextStateTime::MVector{T,Float64},nextInputTime :: MVector{T,Float64} )where{T}   
      
        # which is faster? finding the minimum or implementing a priority queue
        minStateTime=Inf
        minState_index=0  # what if all nextstateTime= Inf ...especially at begining????? min_index stays 0!!!
        minInputTime=Inf
        minInput_index=0
        ST_STATE=1
        ST_INPUT=2
        for i=1:T
            if nextStateTime[i]<minStateTime
                minStateTime=nextStateTime[i]
                minState_index=i
            end
        end
        for i=1:T
            if nextInputTime[i] < minInputTime
                minInputTime=nextInputTime[i]
                minInput_index=i
            end
        end
    
    
    
    #=     if minState_index==0 
            println(" static system! all derivatives are null!")
            return (1,minTime) # later throw exception or maybe draw horizontal lines at initial conditions
        end  =#
        if minInputTime<minStateTime
           # println("an event N",minEvent_index, "about to occur! at time= ",minEventTime)
            return (minInput_index,minInputTime,:ST_INPUT)
        else
            return (minState_index,minStateTime,:ST_STATE)
        end
    
    
      end  
    

 end#end module
 #-------------------------user space-----------------------------------
 using BenchmarkTools
 using Plots;gr()
 function test()
  Debug=false
 odeprob=module1.@ode begin
    u=[1.0,2.4]
    d=[1.0,0.5]
     dq[1]=q[2]#+2.0q[2]*q[1]   
     dq[2]=-q[1]-q[2]
    # du[2]=u[2]-u[1]+cos(t)
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
 sol=module1.integrate(odeprob)
 if !Debug
  temp1=[]
  temp2=[]
  for i=1:length(sol[2][1])
      push!(temp1,sol[2][1][i].coeffs[1])
      push!(temp2,sol[2][2][i].coeffs[1])
  end
  display(plot!(sol[1],temp1))
  display(plot!(sol[1],temp2)) 
  end
  
 

 end

 #@btime test()
 test()
#=  using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2]#*cos(t) #*30*exp(t)
        du[2] = -u[1] - u[2]#+sqrt(t) # +1#+10-t*t+t +cos(t)
       # du[2]=(1 - u[1]^2) * u[2] - u[1] 
      # du[2]=1/(t+1) + sin(t)*sqrt(t)
     # du[2]=sqrt(t)
    
    end
    tspan = (0.0,5)
    u0 = [1.0,2.4]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end
odeDiffEquPackage() =#
 

 
  #@btime modqss.qss_Integrate(initCond,jacobian,order)
  
  
  println("done") 
  readline()
