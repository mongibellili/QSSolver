#using TimerOutputs
#using InteractiveUtils

function QSS_integrate(::Val{2}, s::QSS_data, settings::SimSetting, odep::NLODEProblem)
 # Debug = false
  # Debug=true
  allowPrint=-1e12# for printing after an event
 # reset_timer!()
 #*********************************settings*****************************************
 ft = settings.finalTime
 initTime = settings.initialTime
 relQ = settings.dQrel
 absQ = settings.dQmin
 solver = settings.solver
 savetimeincrement=settings.savetimeincrement
   #*********************************qss method data*****************************************
   quantum = s.quantum
   order=s.order
   nextStateTime = s.nextStateTime
   nextEventTime = s.nextEventTime
   nextInputTime = s.nextInputTime
   tx = s.tx
   tq = s.tq
  #=  x = s.x
   q = s.q =#
   x=Taylor1{Float64}[]
   q=Taylor1{Float64}[]
  #*********************************problem info*****************************************
  initConditions = odep.initConditions
  d = odep.discreteVars
  #----------to compute derivatives
  jacobian = odep.jacobian
  discJac = odep.discreteJacobian
  #inputVars=odep.inputVars
  #----------to compute ZC expressions
  zc_jac = odep.ZC_jacobian
  ZC_discJac = odep.ZC_jacDiscrete
  #ZC_input=odep.ZCinputVars
  #-----------to execute event Dependencys
  evDep = odep.eventDependencies
#=   println("jacobian= ",jacobian)
  println("discJac= ",discJac)
  println("zc_jac= ",zc_jac)
  println("ZC_discJac= ",ZC_discJac)
  println("evDep= ",evDep) =#
 #********************************helper values*******************************

  #numCont = computeStates(initConditions)
  numDiscr = length(d)
  states = length(initConditions)
  numberZC=size(zc_jac, 1)
  numEvents=numberZC*2
  savetime = savetimeincrement
  oldsignValue = MMatrix{numberZC,2}(zeros(numberZC*2))  #usedto track if zc changed sign; each zc has a value and a sign 



 #*******************************create dependencies**************************$$
  SD = createDependencyMatrix(jacobian)
  dD =createDependencyMatrix(discJac) # temp dependency to be used to determine HD1 and HZ1 HD=Hd-dD Union Hs-sD
  SZ =createDependencyMatrix(zc_jac) 
  dZ =createDependencyMatrix(ZC_discJac) # temp dependency to be used to determine HD2 and HZ2
  HZ1HD1=createDependencyToEventsDiscr(dD,dZ,evDep) 
  HZ2HD2=createDependencyToEventsCont(SD,SZ,evDep) 
  HZ=unionDependency(HZ1HD1[1],HZ2HD2[1])
  HD=unionDependency(HZ1HD1[2],HZ2HD2[2])
  

#=   println("SD= ",SD)
  println("SZ= ",SZ)
  println("HZ= ",HZ)
  println("HD= ",HD) =#








  f = Vector{Function}()
  for i = 1:length(odep.eqs)# later change to states
    # push!(arr,eval(eqs[i])) #f_1
   # println("odep.eqs[i].args[2])= ",odep.eqs[i].args[2])
    push!(f, @RuntimeGeneratedFunction(odep.eqs[i].args[2])) #args[2] cuz there is extra stuff (line number in arg1)
  end

  zcf = Vector{Function}()
  for i = 1:length(odep.zceqs)# later change to numberZC
    # push!(arr,eval(eqs[i])) #f_1
    push!(zcf, @RuntimeGeneratedFunction(odep.zceqs[i].args[2])) #args[2] cuz there is extra stuff
  end
#display(zcf);println()
eventf = Vector{Function}()

for i = 1:length(odep.eventEqus)# later change to numberZC
  # push!(arr,eval(eqs[i])) #f_1
 # println("odep.eventEqus[i].args[2])= ", odep.eventEqus[i])
  push!(eventf, @RuntimeGeneratedFunction(odep.eventEqus[i])) #args[2] not needed because the way the expr was constructed in the qssnlproblem is different
end


 
  
 
#*********************************  initialize          *****************************************
#=   x = Taylor1{Float64}[]
  q = Taylor1{Float64}[] =#
  t = Taylor1(order)
  cacheT=Taylor1(order)
  #t=Taylor1([1e-120,1,0],order-1) # this is bad cuz 1/t explodes
  #  initTime=5.562684646268003e-308 ########**********sqrt and xreal exponents complain about expansion around zero....this is bad cuz solution explodes (tested)
  #  initTime= 1e-152 # only a not too small number works
  
  interpolationOrder = 1 # 1 means savedVars = x + t derX
  arr = []
  for i = 1:states
    push!(arr, [])
  end
  savedVars = SVector{states,Array{Taylor1{Float64}}}(tuple(arr...))
  savedTimes = Array{Float64}([initTime])



  #*************************************initialize************************************
  for i = 1:states
    push!(x, Taylor1(zeros(order + 1), order) + initConditions[i])
    push!(q, Taylor1(zeros(order), order - 1))
    push!(savedVars[i], Taylor1(zeros(interpolationOrder + 1), interpolationOrder))
  end
#=   if Debug
    println("initialsavedVars  = ", savedVars)
  end =#
  for k = 1:order # compute initial derivatives for x and q (similar to a recursive way )
    for i = 1:states
      q[i].coeffs[k] = x[i].coeffs[k]  # q computed from x and it is going to be used in the next x
    end
    for i = 1:states
      x[i].coeffs[k+1] = ((differentiate(f[i](q,d, t + initTime), k - 1)).coeffs[1]) / factorial(k) # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
    end
  end
#=   if Debug
    println("initialX with derivatives= ", x)
    println("initialq with derivatives= ", q)
  end =#
#println("f= ",f)
  for i = 1:states
    tx[i] = initTime
    tq[i] = initTime
    for j = 1:interpolationOrder
      savedVars[i][1].coeffs[j] = x[i].coeffs[j]
    end
    quantum[i] = relQ * abs(x[i].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)

    if quantum[i] < absQ
      quantum[i] = absQ
    end
    computeNextTime(Val(2), i, initTime, nextStateTime, x, quantum)
    computeNextInputTime(Val(2), i, initTime, 0.0,f[i](q,d,initTime + t), nextInputTime, x,  quantum)
  end
  for i=1:numberZC
    #println("event i= ", i)
    #= println(q[1])
    println("function return taylor= ",zcf[i](q,d,initTime + t)) =#
    output=zcf[i](x,d,initTime + t).coeffs[1] #test this evaluation
   # println("output= ",output)
    oldsignValue[i,2]=output #value
    oldsignValue[i,1]=sign(output) #sign modify 
    
  end
 # counter=MVector{1,Int}(6) # delete later ...used for limited println 
  for i=1:numberZC
    computeNextEventTime(i,zcf[i](x,d,initTime + t),oldsignValue,initTime,  nextEventTime, quantum)
  end



#=   if Debug
    println("initial quantum= ", quantum)
    println("intial nextstateTime= ", nextStateTime)
    println("intial nextInputTime= ", nextInputTime)
  end =#
  simt = initTime
  count = 0
  ###################################################################################################################################################################
  ####################################################################################################################################################################
  #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
    ###################################################################################################################################################################
  ####################################################################################################################################################################

  while simt < ft #&& count < 6
   # count += 1

    sch = updateScheduler(nextStateTime,nextEventTime, nextInputTime)
   # sch = updateScheduler(nextStateTime,nextEventTime)
    simt = sch[2]
    index = sch[1]

    #  simt, index = findmin(nextStateTime)
#=     if Debug
      println("loop$count : nextStateTime rightafter findmin = ", nextStateTime)
      println("loop$count : nextInputTime rightafter findmin = ", nextInputTime)
      println("loop$count :  ", sch[3])
    end
    if Debug
      println("simtime= ", simt, " index=  ", index)
    end =#
    # println(x[index](t))
   #=  if allowPrint>-2
      println("index=========",index)
      allowPrint-=1
    end =#
  
    elapsed = simt - tx[index]
    #  if sch[3] ==:ST_STATE


    #shift taylor
  #=    for k = 1:order #evaluate x at new time, and derivatives...remember coeffs store the over facto also              
      x[index].coeffs[k] = (((differentiate(x[index], k - 1)(t + elapsed))).coeffs[1]) / factorial(k - 1)
    end =#
    if sch[3] == :ST_STATE
   # x[index]=x[index](t+elapsed)
   integrateState(Val(2),x[index],cacheT,elapsed)
    #if Debug
     # println("loop $count :x after evaluatex = ", x)
   # end
    quantum[index] = relQ * abs(x[index].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
    if quantum[index] < absQ
      quantum[index] = absQ
    end
#=     if Debug
      println("loop $count :quantum = ", quantum)
    end =#
    tx[index] = simt
      for k = 1:order
      q[index].coeffs[k] = x[index].coeffs[k] #updateQ
    end
#=     if Debug
      println("loop $count :q after updateQ = ", q)
    end =#
    tq[index] = simt #tq needed for higher orders
    ##################################################state###################################################
    
   #=    if allowPrint>0
        println("==========   ============   ==========from integrate---at begining of state change:===========   ============   ==========")
        println("at time $simt state change  index= ",index)
        println("nextEventTime= ",nextEventTime)
        println("nextInputTime= ",nextInputTime)
        println("nextstateTime= ",nextStateTime)
        println("at begining x= ",x)
        #allowPrint-=1
      end =#
     # println("at time $simt state change  index= ",index)
       computeNextTime(Val(2), index, simt, nextStateTime, x, quantum) #
      for i = 1:length(SD[index])
        j = SD[index][i] # this line before the next line or vice versa gave the same bench results
        if j != 0
       #    println("dependency loop $j= ")              
          elapsed = simt - tx[j]
          # @timeit "state-integrate state" 
          if elapsed > 0
            #@timeit "integrate state" 
            #x[j].coeffs[1] = x[j](t + elapsed).coeffs[1] # "evaluate" x at new time only...derivatives get updated next using computeDerivative()
            x[j].coeffs[1] = x[j](elapsed)
            q[j].coeffs[1] = x[j].coeffs[1]
            tx[j] = simt
            tq[j] = simt

          end
         # computeDerivative(Val(2), j,tq,x,q,d,t, f[j],cacheT,simt)
           computeDerivative1(Val(2), j,tq,x,q,d,t, f[j],cacheT,simt)
         reComputeNextTime(Val(2), j, simt, nextStateTime, x, q, quantum)
        #=  if allowPrint>0
          
          println("end of statechange:inside SD x[$j]= ",x)
          
          
        end =#
        end#end if j!=0
      end#end for

      for i = 1:length(SZ[index])
        j = SZ[index][i] # this line before the next line or vice versa gave the same bench results
         if j != 0             
          #normally and later i should update q (integrate q=q+e derQ  for higher orders)
          computeNextEventTime(j,zcf[j](x,d,simt + t),oldsignValue,simt,  nextEventTime, quantum)
         # println("there is a continous var that influenced zc ")
        end  
       end
       #= if allowPrint>0 
        println("nextEventTime= ",nextEventTime)
        println("nextInputTime= ",nextInputTime)
        println("nextstateTime= ",nextStateTime)
        println("==========  ============  ==========end state change:================  ============  ===========")
        allowPrint-=1 
      end =#
      ##################################################################"input########################################""
    elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself      
     #=  if allowPrint>0
        println("from integrate: input change after the event index= ",index)
        allowPrint-=1
      end =#
     # println("from integrate: input change after the event index= ",index)
      computeNextInputTime(Val(2), index, simt, elapsed,  f[index](q,d,simt + t), nextInputTime, x,  quantum)
     # computeDerivative(Val(2), index,x, f[index](q,d,t),cacheT,simt)
      for i = 1:length(SD[index])
        j = SD[index][i] # this line before the next line or vice versa gave the same bench results
        if j != 0
          # println("dependency loop $j= ")              
          elapsed = simt - tx[j]
          #@timeit "input-integrate state" 
           if elapsed > 0
            #@timeit "integrate state" 
            x[j].coeffs[1] = x[j](t + elapsed).coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
           # q[j].coeffs[1] = x[j].coeffs[1]
            tx[j] = simt
          end
          #computeDerivative(Val(2), j, x, f[j](q,d,t),cacheT,simt)
          computeDerivative1(Val(2), j,tq,x,q,d,t,f[j],cacheT,simt)
          reComputeNextTime(Val(2), j, simt, nextStateTime, x, q, quantum)
        end#end if j!=0
      end#end for


      for i = 1:length(SZ[index])
        j = SZ[index][i] # this line before the next line or vice versa gave the same bench results
         if j != 0             
          #normally and later i should update q (integrate q=q+e derQ  for higher orders)
          computeNextEventTime(j,zcf[j](x,d,simt + t),oldsignValue,simt,  nextEventTime, quantum)
         # println("there is a continous var (changed because input) that influenced zc ")
        end  
       end


   #################################################################event########################################
    else
      #first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign to know ev+ or ev-
    #  println("*********************************************from integrate: event occured***********************************")
      modifiedIndex=0
    #=   println("zcf[index](q,d,simt + t)=   ",zcf[index](q,d,simt + t))
      println("zcf[index](q,d,simt + t)(simt)=   ",zcf[index](q,d,simt + t)(simt))
      println("q1=   ",q[1])
      println("index= ",index) =#
     # if sign(zcf[index](q,d,simt + t)(simt))>0
     if (zcf[index](x,d,simt + t).coeffs[1])>0       # sign is not needed here
        modifiedIndex=2*index-1   # the  event that just occured is at  this index
      else
        modifiedIndex=2*index
      end         
      #=
       since i am checking if zc>0 then eventFunct1 else eventFunction2....it is better to store the whole thing as an event function
        that way i can have more complicated events without worrying about parsing...one thing: is it cheaper to store one giant event in a function
        or too many small events...i think one giant event is cheaper.




      =#
      # println("modified index= ",modifiedIndex)
      eventf[modifiedIndex](q,d,simt + t)
      #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
    #  println("after execute event d= ",d)
      nextEventTime[index]=Inf   #not tested
    #  println("inside integrate:event --x before event= ",x)
      for i = 1:length(HD[modifiedIndex]) # care about dependency to this event only
            j = HD[modifiedIndex][i] # this line before the next line or vice versa gave the same bench results
             if j != 0
        #      println("inside integrate:event --the state var affected by event is x$j")
                  elapsed = simt - tx[j]
                  if elapsed > 0  # if event triggere by change of sign and time=now then elapsed=0
                    #@timeit "integrate state" 
              #      println("elapsed should appear when event triggered as scheduled not by sign change = ",elapsed)
                    x[j].coeffs[1] = x[j](t + elapsed).coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
                    q[j].coeffs[1] = x[j].coeffs[1]
                    tx[j] = simt
                  end
              #    println("before derivative event---------------------------------- d= ",d)
                  #computeDerivative(Val(2), j,x, f[j](q(simt + t),d,simt + t)) #updating q as a vector updates all components q[0] q[1]...
                  #but some q is already updated!!!causes error. if q[1] updated at 0.5 then q(t+0.5) makes another change
                 # computeDerivative(Val(2), j,x, f[j](q,d,t),cacheT,simt)
                  computeDerivative1(Val(2), j,tq,x,q,d,t, f[j],cacheT,simt)
                  reComputeNextTime(Val(2), j, simt, nextStateTime, x, q, quantum)
                #=   println("inside integrate:end HD  after event x= ",x)
                  println("inside integrate:end HD nexttime= ",nextStateTime)
                  println("inside integrate:end HD nexteventtime= ",nextEventTime) =#
                  
              end
            
      end
      for i = 1:length(HZ[modifiedIndex])
            j = HZ[modifiedIndex][i] # this line before the next line or vice versa gave the same bench results
             if j != 0             
              #normally and later i should update q (integrate q=q+e derQ  for higher orders)
              computeNextEventTime(j,zcf[j](x,d,simt + t),oldsignValue,simt,  nextEventTime, quantum)
            #  println("there is a discrete var that influenced zc ")
            end  
          
      end
    #  println("***********************************************end of event*****************************************************")
     allowPrint+=3
    end#end state input event

    
  
      if simt > savetime
        savetime += savetimeincrement
          for k = 1:states
                temp = Taylor1(interpolationOrder)
                for j = 1:interpolationOrder
                  temp.coeffs[j] = x[k].coeffs[j]
                end
                # println("temp= ",temp)
                push!(savedVars[k], temp)
            end
           push!(savedTimes, simt)
          end#end if save
  end#end while
 # print_timer()
  (savedTimes, savedVars)
end#end integrate