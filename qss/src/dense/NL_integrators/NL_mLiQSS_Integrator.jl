 #using TimerOutputs
 #using TimerOutputs

 function mLiQSS_integrate(CommonqssData::CommonQSS_data{O,0},liqssdata::LiQSS_data{O,false},specialLiqssData::SpecialLiqssQSS_data, odep::NLODEProblem{PRTYPE,T,0,0,CS},f::Function,jac::Function,SD::Function,map::Function) where {PRTYPE,O,T,CS}
  #= cacheA=specialLiqssData.cacheA
  direction=specialLiqssData.direction
  qminus= specialLiqssData.qminus
  buddySimul=specialLiqssData.buddySimul =#
  ft = CommonqssData.finalTime;initTime = CommonqssData.initialTime;relQ = CommonqssData.dQrel;absQ = CommonqssData.dQmin;maxErr=CommonqssData.maxErr;
  #= savetimeincrement=CommonqssData.savetimeincrement;savetime = savetimeincrement =#
  quantum = CommonqssData.quantum;nextStateTime = CommonqssData.nextStateTime;nextEventTime = CommonqssData.nextEventTime;nextInputTime = CommonqssData.nextInputTime
  tx = CommonqssData.tx;tq = CommonqssData.tq;x = CommonqssData.x;q = CommonqssData.q;t=CommonqssData.t
   savedVars=CommonqssData.savedVars;
  savedTimes=CommonqssData.savedTimes;integratorCache=CommonqssData.integratorCache;taylorOpsCache=CommonqssData.taylorOpsCache;cacheSize=odep.cacheSize
  
  #a=deepcopy(odep.initJac);
  a=liqssdata.a
  u=liqssdata.u
  #***************************************************************  
  qaux=liqssdata.qaux;olddx=liqssdata.olddx;olddxSpec=liqssdata.olddxSpec

  numSteps = Vector{Int}(undef, T)
   #######################################compute initial values##################################################
   n=1
  for k = 1:O # compute initial derivatives for x and q (similar to a recursive way )
    n=n*k
      for i = 1:T
        q[i].coeffs[k] = x[i].coeffs[k]  # q computed from x and it is going to be used in the next x
      end
      for i = 1:T
         clearCache(taylorOpsCache,Val(CS),Val(O));f(i,q,t ,taylorOpsCache)
        ndifferentiate!(integratorCache,taylorOpsCache[1] , k - 1)
        x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
      end
  end

 
   for i = 1:T
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

    numSteps[i]=0
    #= @timeit "savevars" =# push!(savedVars[i],x[i][0])
     push!(savedTimes[i],0.0)
     quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
    nupdateQ(Val(O)#=  =#,i,x,q,quantum,a,u,qaux,olddx,tx,tq,initTime,ft,nextStateTime) 
  end
  for i = 1:T
     clearCache(taylorOpsCache,Val(CS),Val(O));f(i,q,t,taylorOpsCache);
     computeDerivative(Val(O), x[i], taylorOpsCache[1])#0.0 used to be elapsed...even down below not neeeded anymore
    Liqss_reComputeNextTime(Val(O), i, initTime, nextStateTime, x, q, quantum,a)
    computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)#not complete, currently elapsed=0.1 is temp until fixed
   
  end
  ###################################################################################################################################################################
  ####################################################################################################################################################################
  #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
  ###################################################################################################################################################################
  #################################################################################################################################################################### 
  simt = initTime ;simulStepCount=0;totalSteps=0;prevStepTime=initTime
 
  simul=false
   while simt < ft && totalSteps < 200000000
    sch = updateScheduler(Val(T),nextStateTime,nextEventTime, nextInputTime)
    simt = sch[2]
    #= if  simt>ft  
      saveLast!(Val(T),Val(O),savedVars, savedTimes,saveVarsHelper,ft,prevStepTime, x)
      break   ###################################################break##########################################
    end =#
    index = sch[1]
    totalSteps+=1
    t[0]=simt
    ##########################################state########################################
    if sch[3] == :ST_STATE
        elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
        quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index] 
           
        nupdateQ(Val(O),index,x,q,quantum,a,u,qaux,olddx,tx,tq,simt,ft,nextStateTime) ;tq[index] = simt        
       # Liqss_ComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum)
        #----------------------------------------------------check dependecy cycles---------------------------------------------      
        for j in SD(index) 
            if j!=index && a[index][j]*a[j][index]!=0.0            
              if isCycle_and_simulUpdate(Val(O),index,j,x,q,quantum,a,u,qaux,olddx,#= olddxSpec ,=#tx,tq,simt,ft)
                simulStepCount+=1             
                Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)
                for k in SD(j)
                        elapsedx = simt - tx[k]
                        if elapsedx > 0
                          x[k].coeffs[1] = x[k](elapsedx);tx[k] = simt
                          integrateOldx(Val(O),x[k],olddx[k],elapsedx)
                          elapsedq = simt - tq[k]
                          if elapsedq > 0 integrateState(Val(O-1),q[k],elapsedq); tq[k] = simt end
                        end
                       
                        for b in (jac(j)  )   
                            elapsedq = simt - tq[b]
                            if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
                        end                     
                        clearCache(taylorOpsCache,Val(CS),Val(O));f(k,q,t,taylorOpsCache);computeDerivative(Val(O), x[k], taylorOpsCache[1])
                        Liqss_reComputeNextTime(Val(O), k, simt, nextStateTime, x, q, quantum,a)
                        updateOtherApprox(k,j,x,q,a,qaux,olddx)
                end#end for k depend on j
                                          
                updateLinearApprox(j,x,q,a,qaux,olddx)             
              end#end ifcycle check
        end#end if j!=0
      end#end FOR_cycle check
      #-------------------------------------------------------------------------------------
      #---------------------------------normal liqss: proceed--------------------------------
      #-------------------------------------------------------------------------------------
  
      for j in SD(index)
          #flag[j]=0.0  # all vars that might triggered flag , return to normal       
          elapsedx = simt - tx[j]
          if elapsedx > 0
            x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt
            integrateOldx(Val(O),x[j],olddx[j],elapsedx) # if elapsedx>0 then elapsedq>0 (confirm?)
            elapsedq = simt - tq[j]
            if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt  end
          else
            #olddx[j][1]=x[j][1]
          end

          #= for b = 1:T # elapsed update all other vars that this derj depends upon
            if jac[j][b] != 0   =# 
          for b in (jac(j)  )    
              elapsedq = simt - tq[b]
              if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
          end
          clearCache(taylorOpsCache,Val(CS),Val(O));f(j,q,t,taylorOpsCache)
            computeDerivative(Val(O), x[j], taylorOpsCache[1])
            Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)
            updateOtherApprox(j,index,x,q,a,qaux,olddx)# Later inverstigate oldspec not updated when aij*aji=0 above


      end#end for SD
   
      # if abs(a[index][index])>1e-6  # if index depends on itself update, otherwise leave zero 
      updateLinearApprox(index,x,q,a,qaux,olddx)########||||||||||||||||||||||||||||||||||||liqss|||||||||||||||||||||||||||||||||||||||||
    #  end
  
      ##################################input########################################
    elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
    
      elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
      quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
      for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt 
        for b in jac(index) 
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
        end
      clearCache(taylorOpsCache,Val(CS),Val(O));f(index,q,t,taylorOpsCache)
      computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
      computeDerivative(Val(O), x[index], taylorOpsCache[1])
  
      for j in (SD(index))    
          elapsedx = simt - tx[j];
          if elapsedx > 0 
            x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt 
            quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]   
          end
          elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext                 
          for b in jac(j)  
              elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
           
          end
          clearCache(taylorOpsCache,Val(CS),Val(O));f(j,q,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
          reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end for
      
    #################################################################event########################################
  
    end#end state/input/event
  
   
 
    push!(savedVars[index],x[index][0])
    push!(savedTimes[index],simt)
 
  end#end while

 #=  for i=1:T# throw away empty points
    resize!(savedVars[i],saveVarsHelper[1])
  end
  resize!(savedTimes,saveVarsHelper[1]) =#

 # print_timer()

 #@timeit "createSol" 
 createSol(Val(T),Val(O),savedTimes,savedVars, "mliqss$O",string(odep.prname),absQ,totalSteps,simulStepCount,numSteps,ft)
     # change this to function /constrcutor...remember it is bad to access structs (objects) directly
  end#end integrate
      
      
