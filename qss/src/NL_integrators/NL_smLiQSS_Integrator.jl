 #using TimerOutputs
function smLiQSS_integrate(::Val{O}, s::LiQSS_data{T,Z,O}, odep::NLODEProblem{T,D,Z,Y},f::Function) where {O,T,D,Z,Y}
  #reset_timer!()
    #*********************************settings*****************************************
    #printCounter=[0,0]#vector{Int} fort debugging to be deleted
    ft = s.finalTime
    initTime = s.initialTime
    relQ = s.dQrel
    absQ = s.dQmin
    savetimeincrement=s.savetimeincrement
    #*********************************qss method data*****************************************
    quantum = s.quantum
    nextStateTime = s.nextStateTime
    nextEventTime = s.nextEventTime
    nextInputTime = s.nextInputTime
    tx = s.tx
    tq = s.tq
    x = s.x
    q = s.q
    t=s.t
    savedVars=s.savedVars
    savedTimes=s.savedTimes
    integratorCache=s.integratorCache
    taylorOpsCache=s.taylorOpsCache
    cacheSize=odep.cacheSize
    #*********************************problem info*****************************************
    initConditions = odep.initConditions
    d = odep.discreteVars
    #----------to compute derivatives
    jacobian = odep.jacobian
    jac=changeBasicToInts(jacobian)# change from type nonisbits to int so that access is cheaper down
    a=s.initJac
    u=s.u
    tu=s.tu
    qaux=s.qaux
    olddx=s.olddx

    discJac = odep.discreteJacobian
    #----------to compute ZC expressions
    zc_jac = odep.ZC_jacobian
    ZC_discJac = odep.ZC_jacDiscrete
    #-----------to execute event Dependencys
    evDep = odep.eventDependencies
    #********************************helper values*******************************
    numDiscr = length(d)
    #states = length(initConditions)
    #numberZC=size(zc_jac, 1)
    numEvents=Z*2   #each zero crossing has 1 posEv and 1 negEv
    savetime = savetimeincrement
    oldsignValue = MMatrix{Z,2}(zeros(Z*2))  #usedto track if zc changed sign; each zc has a value and a sign 
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
    zcf = Vector{Function}()
    for i = 1:length(odep.zceqs)# later change to numberZC
      push!(zcf, @RuntimeGeneratedFunction(odep.zceqs[i].args[2])) #args[2] cuz there is extra stuff
    end
    eventf = Vector{Function}()
    for i = 1:length(odep.eventEqus)# later change to numEvents
      push!(eventf, @RuntimeGeneratedFunction(odep.eventEqus[i].args[2])) 
    end
    
    #######################################compute initial values##################################################
    n=1
    for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
      n=n*k
       for i = 1:T
          q[i].coeffs[k] = x[i].coeffs[k]  # q computed from x and it is going to be used in the next x
        end
        for i = 1:T
          clearCache(taylorOpsCache,cacheSize)
           f(i,q,d, t ,taylorOpsCache)
           (ndifferentiate!(integratorCache,taylorOpsCache[1] , k - 1))
          x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
        end
    end
    
    for i = 1:T
      #= for k=1:O# deleting this causes scheduler error
          u[i][k]=x[i][k]-q[i][k-1]*a[i][i] #  later we will investigate inconsistencies of using data stored vs */ factorial!!! ...also do not confuse getindex for taylor...[0] first element and u[i][1]...first element    ########||||||||||||||||||||||||||||||||||||liqss|||||||||||||||||||||||||||||||||||||||||
      end =#

      #= for j=1:T
        if j!=i
          u[i][j][1]=x[i][1]-a[i][i]*q[i][0]-a[i][j]*q[j][0]
        else
          u[i][j][1]=x[i][1]-a[i][i]*q[i][0]
        end
      end =#

      p=1
      for k=1:O
        p=p*k
        m=p/k
        for j=1:T
          if j!=i
            u[i][j][k]=p*x[i][k]-a[i][i]*m*q[i][k-1]-a[i][j]*m*q[j][k-1]
          else
            u[i][j][k]=p*x[i][k]-a[i][i]*m*q[i][k-1]
          end
        end
      end
      savedVars[i][1].coeffs .= x[i].coeffs  #to be changed  1 2 3 ?
      quantum[i] = relQ * abs(x[i].coeffs[1]) 
      if quantum[i] < absQ
        quantum[i] = absQ
      end

      updateQ(Val(O),i,x,q,quantum,a,u,qaux,olddx,tq,tu,initTime,ft,nextStateTime) 
     
        
         clearCache(taylorOpsCache,cacheSize)
         f(i,q,d,t,taylorOpsCache) 
   
        
         computeDerivative(Val(O), x[i], taylorOpsCache[1],integratorCache,0.0)#0.0 used to be elapsed...even down below not neeeded anymore
        
          Liqss_reComputeNextTime(Val(O), i, initTime, nextStateTime, x, q, quantum,a)




      computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)#not complete, currently elapsed=0.1 is temp until fixed
    end
   # @show u
    for i=1:Z
      clearCache(taylorOpsCache,cacheSize)
      output=zcf[i](x,d,t,taylorOpsCache).coeffs[1] #test this evaluation
      oldsignValue[i,2]=output #value
      oldsignValue[i,1]=sign(output) #sign modify 
      computeNextEventTime(i,output,oldsignValue,initTime,  nextEventTime, quantum)#,printCounter)# if Z=0 nexteventtime did not get filled with zeros in the first place
    end
    ###################################################################################################################################################################
    ####################################################################################################################################################################
    #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
    ###################################################################################################################################################################
    ####################################################################################################################################################################
    simt = initTime
    count = 1 # not zero because intial value took 0th position
    len=length(savedTimes)
    printcount=0
    while simt < ft #&& printcount < 20000
      printcount+=1
      sch = updateScheduler(nextStateTime,nextEventTime, nextInputTime)
      simt = sch[2]
     #=  if 10<simt<50.1
        @show simt
        @show printcount
      end =#
      index = sch[1]
      if debug
      @show simt
      @show index
      end
      t[0]=simt
      ##########################################state########################################
      if sch[3] == :ST_STATE
        elapsed = simt - tx[index]
        integrateState(Val(O),x[index],integratorCache,elapsed)
                                                                              if debug  println("x after intgrate ", x[index]) end
        quantum[index] = relQ * abs(x[index].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
        if quantum[index] < absQ
          quantum[index] = absQ
        end
                                                                               if debug  println("quantum$index= ",quantum[index] ) end
        #= mupdateQ(Val(O),index,x,q,quantum,a,u,qaux,olddx,tx,tq,tu,simt,ft) ########
        computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) # =#
        updateQ(Val(O),index,x,q,quantum,a,u,qaux,olddx,tq,tu,simt,ft,nextStateTime) 
        tx[index] = simt
        tq[index] = simt #done inside : either way is fine
                                                                                    if debug println("tq[$index]= ",tq[index]) end
       #----------------------------------------------------check dependecy cycles---------------------------------------------    
       for l = 1:length(SD[index])
        j = SD[index][l] 
        if j != 0 && j!=index && a[index][j]*a[j][index]!=0           
          elapsed = simt - tx[j]
          if isCycle_and_simulUpdate(Val(O),index,j,x,q,quantum,a,u,qaux,olddx,tx,tq,tu,simt,ft)
                    #if 10.3<simt<12 
             ###         println("-------after if iscycle nextstate was = ",nextStateTime)
                   # end
                    Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)
                   # if 10.3<simt<12 
                  
              ###        println("-                 nextstatetime became ",nextStateTime)
                   # end
                  #  Liqss_reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum,a)
               for l = 1:length(SD[j])
                      k = SD[j][l] 
                      if k != 0  
                        
                        elapsedx = simt - tx[k]
                        if elapsedx > 0
                          #"evaluate" x only at new time ...derivatives get updated next using computeDerivativ()
                          x[k].coeffs[1] = x[k](elapsedx)
                         
                          
                                                                 #  if debug println("x$j elapse updated")  end 
                          tx[k] = simt
                        
                        end
                        elapsedq = simt - tq[k]
                        if elapsedq > 0
                          #"evaluate" x only at new time ...derivatives get updated next using computeDerivativ()
                          
                          integrateState(Val(O-1),q[k],integratorCache,elapsedq)
                         # q[j].coeffs[1] = q[j](elapsedq) # ouch ! this bit me for a day: q needs to be updated here for recomputeNext, not just for the derivatives!!!
                                                                  # if debug println("x$j elapse updated")  end 
                          
                         tq[k] = simt
                        end
        
                        for b = 1:T # elapsed update all other vars that this derj depends upon.needed for when sys has 3 or more vars.
                          #sj = jac[k][b] 
                          if jac[k][b] != 0  
                                                               #  if debug  @show sj   end      
                            elapsedq = simt - tq[b]
                            if elapsedq>0
                             # q[b].coeffs[1] = q[b](elapsed) ## 
                              integrateState(Val(O-1),q[b],integratorCache,elapsedq)
                              tq[b]=simt
                                                                 # if debug println("q$b elapse updated under sj") end
                            end
                          end
                        end



                       #=  elapsed = simt - tx[k]
                        if elapsed>0
                         x[k].coeffs[1] = x[k](elapsed) #
                         #olddx[k][1]=x[k][1]+2*x[k][2]*elapsed #olddx always has to be elapsed updated
                         updateOlddx(Val(O), k,x,olddx,elapsed)
                         tx[k]=simt
                        end
                                                                                            if debug
                                                                                            println("k depends on j; k=",k)
                                                                                            @show x[k][0],q[k][0],olddx[k][1]
                                                                                            end
                          for b = 1:length(SD[k]) # elapsed update all other vars that this derK depends upon
                            s = SD[k][b] 
                            if s != 0           
                              elapsed = simt - tq[s]
                              if elapsed>0
                                #q[s].coeffs[1] = q[s](elapsed) #
                                integrateState(Val(O-1),q[s],integratorCache,elapsed)
                                tq[s]=simt
                              end
                            end
                          end =#
                        clearCache(taylorOpsCache,cacheSize)
                        f(k,q,d,t,taylorOpsCache)
                        computeDerivative(Val(O), x[k], taylorOpsCache[1],integratorCache,elapsed)
                      # @timeit "state-recompute" 
                    #  if 10.3<simt<12 
                        
                     #   @show nextStateTime
                     # end
                        Liqss_reComputeNextTime(Val(O), k, simt, nextStateTime, x, q, quantum,a)
                  ###      println("-*-*-*-*-*inside k dependcy nextStateTime became= ",nextStateTime)
                        
                        
                                                                      if debug 
                                                                        @show x[k][0],q[k][0],x[k][1]
                                                                        println("begining of update other after cycle detected: akj= ",a[k][j])
                                                                      end
                        updateOtherApprox(Val(O),k,j,x,q,a,u,qaux,olddx,tu,simt)
                                                                       if debug println("after akj update: a$k$j= ",a[k][j]) end    
                      end#end if k!=0
                end#end for k depend on j
                for l = 1:length(SZ[j])
                      k = SZ[index][l] 
                      if k != 0             
                        #normally and later i should update q (integrate q=q+e derQ  for higher orders)
                        clearCache(taylorOpsCache,cacheSize)
                        computeNextEventTime(k,zcf[k](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
                      end  #end if j!=0
                 end#end for SZ
                                                                             if debug  println("begining of update other after cycle detected: ajj= ",a[j][j]) end
                 updateLinearApprox(Val(O),j,x,q,a,u,qaux,olddx,tu,simt)
                                                                             if debug println("after ajj update: a$j$j= ",a[j][j]) end
                 
          end#end ifcycle check
         # tx[j] = simt #  for sys with 2 vars: i won't be updated because of elapsed and j wont be updated
         # tq[j] = simt
        end
      end#end FOR_cycle check
       #---------------------------------normal liqss: proceed--------------------------------
                                                                       if debug println("**********normal dependency**************") end
       for l = 1:length(SD[index])
          j = SD[index][l] 
                                                                          if debug println("j depends on i; j=",j) end
          if j != 0           
            elapsedx = simt - tx[j]
            if elapsedx > 0
              #"evaluate" x only at new time ...derivatives get updated next using computeDerivativ()
              x[j].coeffs[1] = x[j](elapsedx)
             
              
                                                       if debug println("x$j elapse updated")  end 
              tx[j] = simt
            
            end
            elapsedq = simt - tq[j]
            if elapsedq > 0
              #"evaluate" x only at new time ...derivatives get updated next using computeDerivativ()
              
              integrateState(Val(O-1),q[j],integratorCache,elapsedq)
             # q[j].coeffs[1] = q[j](elapsedq) # ouch ! this bit me for a day: q needs to be updated here for recomputeNext, not just for the derivatives!!!
                                                       if debug println("x$j elapse updated")  end 
              
             tq[j] = simt
            end

            for b = 1:T # elapsed update all other vars that this derj depends upon.needed for when sys has 3 or more vars.
              #sj = jacobian[j][b] 
              if jac[j][b] != 0  
                                                     if debug  @show sj   end      
                elapsedq = simt - tq[b]
                if elapsedq>0
                 # q[b].coeffs[1] = q[b](elapsed) ## 
                  integrateState(Val(O-1),q[b],integratorCache,elapsedq)
                  tq[b]=simt
                                                      if debug println("q$b elapse updated under sj") end
                end
              end
            end
           clearCache(taylorOpsCache,cacheSize)
            f(j,q,d,t,taylorOpsCache)
         #   @timeit "comp der"
             computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
                                                                                              if debug println(" x[$j][1]= ", x[j][1]) end
           #=  if 10.3<simt<12 
              println("in normal dependcy nextstate was ",nextStateTime)
              
            end =#
             Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)
             #= if 10.3<simt<12 
              println("in normal dependcy nextstate became ",nextStateTime)
              @show simt
              if nextStateTime[1]>29 || nextStateTime[2]>29  || simt==10.372457046292588

                @show x
                @show q
                @show quantum
                @show a
              end
            end =#
             updateOtherApprox(Val(O),j,index,x,q,a,u,qaux,olddx,tu,simt)
                                                                                         if debug  println("a$j$index after updateoher= ",a[j][index])end
          end#end if j!=0
        end#end for SD
        for l = 1:length(SZ[index])
          j = SZ[index][l] 
          if j != 0             
            #normally and later i should update q (integrate q=q+e derQ  for higher orders)
            clearCache(taylorOpsCache,cacheSize)
            computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
          end  #end if j!=0
        end#end for SZ
        updateLinearApprox(Val(O),index,x,q,a,u,qaux,olddx,tu,simt)########||||||||||||||||||||||||||||||||||||liqss|||||||||||||||||||||||||||||||||||||||||
       
                                                                                          if debug
                                                                                            println("after update aii: a$index$index =",a[index][index]) 
                                                                                            println("tx[$index]= ",tx[index])
                                                                                            println("-----------------------extra verification -------------------------")
                                                                                            var2index=3-index
                                                                                            dxi=a[index][index]*q[index][0]+a[index][var2index]*q[var2index][0]+u[index][var2index][1]
                                                                                            @show dxi,x[index][1]
                                                                                            dxj=a[var2index][var2index]*q[var2index][0]+a[var2index][index]*q[index][0]+u[var2index][index][1]
                                                                                            @show dxj,x[var2index][1]
                                                                                            println("-----------------------end of loop-------------------------")
                                                                                          end
       # tx[index] = simt
        ##################################input########################################
      elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself  

        elapsed = simt - tx[index]    
        clearCache(taylorOpsCache,cacheSize)   
        f(index,q,d,t,taylorOpsCache)
         computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
        for i = 1:length(SD[index])
          j = SD[index][i] 
          if j != 0             
            elapsed = simt - tx[j]
            if elapsed > 0
    
              x[j].coeffs[1] = x[j](elapsed)#.coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivati()
              tx[j] = simt
            end
            quantum[j] = relQ * abs(x[j].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
            if quantum[j] < absQ
              quantum[j] = absQ
            end
            clearCache(taylorOpsCache,cacheSize)
            f(j,q,d,t,taylorOpsCache)
            computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
           # computeDerivative(Val(O), x[j], taylorOpsCache[1])
            Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)
          end#end if j!=0
        end#end for
        for i = 1:length(SZ[index])
          j = SZ[index][i] 
          if j != 0             
            #normally and later i should update q (integrate q=q+e derQ  for higher orders)
            clearCache(taylorOpsCache,cacheSize)
            
            computeNextEventTime(j,zcf[j](q,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter) 
          end  
         # println("end input:who is resizing?")
        end
      #################################################################event########################################
      else
        #first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 

        modifiedIndex=0
        if (zcf[index](x,d,t,taylorOpsCache).coeffs[1])>0       # sign is not needed here
          modifiedIndex=2*index-1   # the  event that just occured is at  this index
        else
          modifiedIndex=2*index
        end       
        eventf[modifiedIndex](q,d,t,taylorOpsCache)
        #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
        nextEventTime[index]=Inf   #investigate more
        for i = 1:length(HD[modifiedIndex]) # care about dependency to this event only
            j = HD[modifiedIndex][i] # 
            if j != 0
                  elapsed = simt - tx[j]             
                  if elapsed > 0  # if event triggere by change of sign and time=now then elapsed=0
                      # println("elapsed appear?= ",elapsed) #only appear when zc depends on x1 and it affects derx2 && derx1 doesn't depend on x1             
                      x[j].coeffs[1] = x[j](elapsed)#.coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivativ()
                      q[j].coeffs[1] = x[j].coeffs[1]
                      tx[j] = simt
                  end
                    clearCache(taylorOpsCache,cacheSize)
    
                  f(j,q,d,t,taylorOpsCache)
                  computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
                 # computeDerivative(Val(O), x[j], taylorOpsCache[1])
                  Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)  
            end     
        end
        for i = 1:length(HZ[modifiedIndex])
              j = HZ[modifiedIndex][i] # this line before the next line or vice versa gave the same bench results
                if j != 0             
                #normally and later i should update q (integrate q=q+e derQ  for higher orders)
                clearCache(taylorOpsCache,cacheSize)           
                computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
              end        
        end
      end#end state/input/event

     #=  @show x
      @show simt =#
      if simt > savetime
          count += 1
          savetime += savetimeincrement #next savetime
          if len<count
            len=count*2
            for i=1:T
              resize!(savedVars[i],len)
              for z=count:len
              savedVars[i][z]=Taylor0(zeros(O+1),O) # without this, the new ones are undefined
              end
            end
            resize!(savedTimes,len)
          end
          for k = 1:T
              
              savedVars[k][count].coeffs .=x[k].coeffs 
          end
         savedTimes[count]=simt
      end#end if save
    end#end while
    for i=1:T# throw away empty points
      resize!(savedVars[i],count)
    end
    #print_timer()
    @show printcount
    resize!(savedTimes,count)
    Sol(savedTimes, savedVars,"mliqss$O")
    end#end integrate
    
    