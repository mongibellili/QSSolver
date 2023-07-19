#using TimerOutputs
function QSS_integrate(CommonqssData::CommonQSS_data{O,T,Z}, odep::NLODEProblem{PRTYPE,T,Z,Y},f::Function,jac::Function,SD::Function,map::Function) where {PRTYPE,O,T,Z,Y}
  
  ft = CommonqssData.finalTime;initTime = CommonqssData.initialTime;relQ = CommonqssData.dQrel;absQ = CommonqssData.dQmin;maxErr=CommonqssData.maxErr;
  savetimeincrement=CommonqssData.savetimeincrement;savetime = savetimeincrement
  quantum = CommonqssData.quantum;nextStateTime = CommonqssData.nextStateTime;nextEventTime = CommonqssData.nextEventTime;nextInputTime = CommonqssData.nextInputTime
  tx = CommonqssData.tx;tq = CommonqssData.tq;x = CommonqssData.x;q = CommonqssData.q;t=CommonqssData.t
   savedVars=CommonqssData.savedVars;
  savedTimes=CommonqssData.savedTimes;integratorCache=CommonqssData.integratorCache;taylorOpsCache=CommonqssData.taylorOpsCache;cacheSize=odep.cacheSize
  prevStepVal = specialLiqssData.prevStepVal
  #*********************************problem info*****************************************
  d = odep.discreteVars
  

  zc_SimpleJac=odep.ZCjac

  HZ=odep.HZ
  HD=odep.HD
  SZ=odep.SZ
 
  evDep = odep.eventDependencies
  #********************************helper values*******************************  

  oldsignValue = MMatrix{Z,2}(zeros(Z*2))  #usedto track if zc changed sign; each zc has a value and a sign 

#######################################compute initial values##################################################
n=1
for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
  n=n*k
   for i = 1:T q[i].coeffs[k] = x[i].coeffs[k] end # q computed from x and it is going to be used in the next x
   for i = 1:T
      clearCache(taylorOpsCache,Val(CS),Val(O));f(i,-1,-1,q,d, t ,taylorOpsCache)
      ndifferentiate!(integratorCache,taylorOpsCache[1] , k - 1)
      x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
    end
end

for i = 1:T
  initSavedVars!(savedVars[i],x[i])
  quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
  computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
  initSmallAdvance=0.1
  #t[0]=initTime#initSmallAdvance
  clearCache(taylorOpsCache,Val(CS),Val(O));
  #@timeit "f" 
  f(i,-1,-1,q,d,t,taylorOpsCache);#@show taylorOpsCache
  computeNextInputTime(Val(O), i, initTime, initTime,taylorOpsCache[1] , nextInputTime, x,  quantum)
  assignXPrevStepVals(Val(O),prevStepVal,x,i)
end

#@show nextStateTime,nextInputTime
for i=1:Z
  #= clearCache(taylorOpsCache,Val(CS),Val(O));output=zcf[i](x,d,t,taylorOpsCache).coeffs[1]  =#
  clearCache(taylorOpsCache,Val(CS),Val(O));
  #@timeit "zcf" 
  f(-1,i,-1,x,d,t,taylorOpsCache)        
                 
  oldsignValue[i,2]=taylorOpsCache[1][0] #value
  oldsignValue[i,1]=sign(taylorOpsCache[1][0]) #sign modify 

  computeNextEventTime(i,taylorOpsCache[1],oldsignValue,initTime,  nextEventTime, quantum)

end

###################################################################################################################################################################
####################################################################################################################################################################
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
###################################################################################################################################################################
####################################################################################################################################################################
simt = initTime ;totalSteps=0;prevStepTime=initTime
  
while simt < ft && totalSteps < 5000000
  sch = updateScheduler(nextStateTime,nextEventTime, nextInputTime)
  simt = sch[2]
 # @timeit "saveLast" 
   if  simt>ft  
    saveLast!(Val(T),Val(O),savedVars, savedTimes,saveVarsHelper,ft,prevStepTime,integratorCache, x)
    break   ###################################################break##########################################
  end
  index = sch[1]
  totalSteps+=1
  t[0]=simt
  ##########################################state######################################## 
  if sch[3] == :ST_STATE
  

  
    elapsed = simt - tx[index];integrateState(Val(O),x[index],integratorCache,elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt    
    computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) #
   # if simt>=3.09e-5
 
  #  end
  for j in (SD(index))
    elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt end
    # quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]         
    elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext        
    for b in (jac(j)  )    
      elapsedq = simt - tq[b]
      if elapsedq>0
        integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt
      end
     end
        clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
      
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
          
    end#end for SD
    for j in (SZ[index])
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
            elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],integratorCache,elapsedx);tx[b]=simt end
            #elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
        end            
        clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,x,d,t,taylorOpsCache)        
        computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
    end#end for SZ
    ##################################input########################################
  elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
    
    elapsed = simt - tx[index];integrateState(Val(O),x[index],integratorCache,elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt 
      for b in jac(index) 
        elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
      end
    clearCache(taylorOpsCache,Val(CS),Val(O));f(index,-1,-1,q,d,t,taylorOpsCache)
    computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
    computeDerivative(Val(O), x[index], taylorOpsCache[1],integratorCache,elapsed)

    for j in(SD(index))  
      elapsedx = simt - tx[j];
      if elapsedx > 0 
        x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt 
        quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]   
      end
      elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext                 
      # elapsed update all other vars that this derj depends upon.
        for b in jac(j) 
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
        end
      
        clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
    end#end for
    for j in (SZ[index])
      
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
             
            elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],integratorCache,elapsedx);tx[b]=simt end
           # elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
       
        end              
       #=  clearCache(taylorOpsCache,Val(CS),Val(O))#normally and later i should update x,q (integrate q=q+e derQ  for higher orders)
        computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer)  =#
        clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,x,d,t,taylorOpsCache)        
         computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
     
    end
  #################################################################event########################################
  else
          printcounter=5
         #=  println("x at start of event")
          @show x
          @show simt 
          @show quantum =#
      
          for b in zc_SimpleJac[index] # elapsed update all other vars that this zc depends upon.
              
              elapsedx = simt - tx[b];
              if elapsedx>0 
                integrateState(Val(O),x[b],integratorCache,elapsedx);
               # println("should update var1 b=",b)
                tx[b]=simt 
              end
             elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
           
          end    
          modifiedIndex=0#first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 
         
          clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,index,-1,x,d,t,taylorOpsCache)    # run ZCF-------- 
         #=  println(" just after event")
         
          @show t 
          @show nextEventTime    =#
                
          if (taylorOpsCache[1][0])>=-1e-10       # sign is not needed here
            modifiedIndex=2*index-1   # the  event that just occured is at  this index
            #println("event1")
          else
            modifiedIndex=2*index
          end  
        # println("x before event= ",x)  
           
         for b in evDep[modifiedIndex].evContRHS
                
                  elapsedx = simt - tx[b];
                  if elapsedx>0 
                    integrateState(Val(O),x[b],integratorCache,elapsedx);
                  #  println("var2 here b = ",b)
                    
                    tx[b]=simt 
                   end
               
           end
           
         ###### eventf[modifiedIndex](x,d,t,taylorOpsCache) #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
        # @show x,modifiedIndex
         f(-1,-1,modifiedIndex,x,d,t,taylorOpsCache)# execute event--------
        
          for i in evDep[modifiedIndex].evCont
            #------------event influences a Continete var: x already updated in event...here update quantum and q and computenext
            
                #quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
                q[i][0]=x[i][0];tx[i] = simt;tq[i] = simt # for liqss updateQ?
              #  computeNextTime(Val(O), i, simt, nextStateTime, x, quantum) 
           
          end
          nextEventTime[index]=Inf   #investigate more 
          for j in (HD[modifiedIndex]) # care about dependency to this event only
                 
              elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt;#= @show j,x[j] =# end
              elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],integratorCache,elapsedq);tq[j] = simt;#= @show q[j] =#  end#q needs to be updated here for recomputeNext                 
              for b = 1:T # elapsed update all other vars that this derj depends upon.
                if b in jac(j)   
                  elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt;#= @show q[b] =# end
                end
              end
              clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsedx)
              reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
              #@show j,x
           
          end
          for j in (HZ[modifiedIndex])
                 
                    for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
                          
                        elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],integratorCache,elapsedx);tx[b]=simt end
                       #elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],integratorCache,elapsedq);tq[b]=simt end
                    
                    end            
                  #= clearCache(taylorOpsCache,Val(CS),Val(O)) #normally and later i should update q (integrate q=q+e derQ  for higher orders)          
                  computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer) =#
                  
                  clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,x,d,t,taylorOpsCache)        
                 computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
               
              # if 0.4>simt > 0.31  println("$index $j nexteventtime from HZ= ",nextEventTime)   end   
          end
         #=  println("x end of step event")
         @show x 
         @show q =# 
  end#end state/input/event
  if simt > savetime || sch[3] ==:ST_EVENT
    save!(Val(O),savedVars , savedTimes , saveVarsHelper,prevStepTime ,simt,tx ,tq , integratorCache,x , q,prevStepVal)
    savetime += savetimeincrement #next savetime 
  else#end if save
    for k = 1:T  
      
        elapsed = simt - tx[k];integrateState(Val(O),x[k],integratorCache,elapsed);tx[k] = simt #in case this point did not get updated.  
        elapsedq = simt - tq[k];integrateState(Val(O-1),q[k],integratorCache,elapsedq);tq[k]=simt        
      
      assignXPrevStepVals(Val(O),prevStepVal,x,k)
    end
  end
  prevStepTime=simt
end#end while

for i=1:T# throw away empty points
  resize!(savedVars[i],saveVarsHelper[1])
end
resize!(savedTimes,saveVarsHelper[1])

createSol(Val(T),Val(O),savedTimes,savedVars, "qss$O",string(nameof(f)),absQ,totalSteps,0)#0 I track simulSteps 
end#end integrate