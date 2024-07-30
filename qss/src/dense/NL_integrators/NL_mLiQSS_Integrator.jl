
function integrate(CommonqssData::CommonQSS_data{SN,O,0},liqssdata::LiQSS_data{Detection},specialLiqssData::SpecialLiqssQSS_data, odep::NLODEProblem{PRTYPE,T,0,0,CS},f::Function,jac::Function,SD::Function,exacteA::Function) where {PRTYPE,SN,CS,Detection,O,T}
  cacheA=specialLiqssData.cacheA
  direction=specialLiqssData.direction
  #qminus= specialLiqssData.qminus
  buddySimul=specialLiqssData.buddySimul
  ft = CommonqssData.finalTime;initTime = CommonqssData.initialTime;relQ = CommonqssData.dQrel;absQ = CommonqssData.dQmin;maxErr=CommonqssData.maxErr;
  #setprecision(BigFloat,64)
  #savetimeincrement=CommonqssData.savetimeincrement;savetime = savetimeincrement
  quantum = CommonqssData.quantum;nextStateTime = CommonqssData.nextStateTime;nextEventTime = CommonqssData.nextEventTime;
  nextInputTime = CommonqssData.nextInputTime
  tx = CommonqssData.tx;tq = CommonqssData.tq;x = CommonqssData.x;q = CommonqssData.q;t=CommonqssData.t
  savedVars=CommonqssData.savedVars;
  savedTimes=CommonqssData.savedTimes;integratorCache=CommonqssData.integratorCache;taylorOpsCache=CommonqssData.taylorOpsCache;#cacheSize=odep.cacheSize
 
  #prevStepVal = specialLiqssData.prevStepVal
  #a=deepcopy(odep.initJac);
  #a=liqssdata.a
  #u=liqssdata.u;#tu=liqssdata.tu
  #***************************************************************  
  qaux=liqssdata.qaux;dxaux=liqssdata.dxaux#= olddx=liqssdata.olddx; ; olddxSpec=liqssdata.olddxSpec =#

  numSteps = Vector{Int}(undef, T)

 #=  temporaryhelper = Vector{Int}(undef, 1)
  temporaryhelper[1]=0 =#
#=   savedVarsQ = Vector{Vector{Float64}}(undef, T)  
  savedVarsQ[1]=Vector{Float64}() ;savedVarsQ[2]=Vector{Float64}()     =#
  acceptedi=Vector{Vector{Float64}}(undef,3);acceptedj=Vector{Vector{Float64}}(undef,3);
for i =1:3
  acceptedi[i]=zeros(2)
  acceptedj[i]=zeros(2)
end
  exacteA(q,cacheA,1,1)  # this 'unnecessary call' 'compiles' the function and it helps remove allocations when used after !!!

  #@show exacteA
   #######################################compute initial values##################################################
  n=1
  for k = 1:1 # compute initial derivatives for x and q (similar to a recursive way )
    n=n*k
      for i = 1:T
        q[i].coeffs[k] = x[i].coeffs[k]  # q computed from x and it is going to be used in the next x
      end
      for i = 1:T
         clearCache(taylorOpsCache,Val(CS),Val(1));f(i,q,t ,taylorOpsCache)
        ndifferentiate!(integratorCache,taylorOpsCache[1] , k - 1)
        x[i].coeffs[k+1] = (integratorCache.coeffs[1]) / n # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
      end
  end

 
   for i = 1:T
    numSteps[i]=0
    #= @timeit "savevars" =# push!(savedVars[i],x[i][0])
    #push!(savedVarsQ[i],q[i][0])
     push!(savedTimes[i],0.0)
     quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
    
    updateQ(Val(O),i,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,initTime,ft,nextStateTime)  
  end
  for i = 1:T
     clearCache(taylorOpsCache,Val(CS),Val(1));f(i,q,t,taylorOpsCache);
     computeDerivative(Val(1), x[i], taylorOpsCache[1])#0.0 used to be elapsed...even down below not neeeded anymore
    Liqss_reComputeNextTime(Val(1), i, initTime, nextStateTime, x, q, quantum)
   
   #prevStepVal[i]=x[i][0]#assignXPrevStepVals(Val(1),prevStepVal,x,i)
  end


  ###################################################################################################################################################################
  ####################################################################################################################################################################
  #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
  ###################################################################################################################################################################
  #################################################################################################################################################################### 
  simt = initTime ;simulStepCount=0;totalSteps=0;
 #totalStepswhenCycles=0
  #simul=false


  while simt < ft && totalSteps < 10000000
    sch = updateScheduler(Val(T),nextStateTime,nextEventTime, nextInputTime) 
    simt = sch[2]
    if simt > ft
      break # 
    end
    index = sch[1]
    numSteps[index]+=1;totalSteps+=1
    t[0]=simt
    ##########################################state########################################
    if sch[3] == :ST_STATE
      xitemp=x[index][0]
        elapsed = simt - tx[index];integrateState(Val(1),x[index],elapsed);tx[index] = simt 
        quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index] 
         #dirI=x[index][0]-savedVars[index][end]  
        
      #=   for b in (jac(index)  )    # update Qb : to be used to calculate exact Aindexb
          elapsedq = simt - tq[b] ;
          if elapsedq>0 integrateState(Val(1-1),q[b],elapsedq);tq[b]=simt end
        end =#
       updateQ(Val(O),index,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,simt,ft,nextStateTime) ;tq[index] = simt   
        #----------------------------------------------------check dependecy cycles---------------------------------------------                
        for j in SD(index)
         #=  for b in (jac(j)  )    # update Qb: to be used to calculate exacte Ajb
            elapsedq = simt - tq[b] ;
            if elapsedq>0  integrateState(Val(1-1),q[b],elapsedq); tq[b]=simt  end
          end =#
          exacteA(q,cacheA,index,j);aij=cacheA[1]
          exacteA(q,cacheA,j,index);aji=cacheA[1]      
          if j!=index && aij*aji!=0.0
              #prvStepValj= savedVars[j][end]#getPrevStepVal(prevStepVal,j) 
              for i =1:3
                acceptedi[i][1]=0.0; acceptedi[i][2]=0.0
                acceptedj[i][1]=0.0; acceptedj[i][2]=0.0
              end 
             
              if isCycle_and_simulUpdate(Val(SN),Val(Detection),acceptedi,acceptedj,aij,aji,index,j,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,simt,ft)
                simulStepCount+=1
               # clearCache(taylorOpsCache,Val(CS),Val(1));f(index,q,t,taylorOpsCache);computeDerivative(Val(1), x[index], taylorOpsCache[1])
              #  clearCache(taylorOpsCache,Val(CS),Val(1));f(j,q,t,taylorOpsCache);computeDerivative(Val(1), x[j], taylorOpsCache[1])
               # Liqss_reComputeNextTime(Val(1), index, simt, nextStateTime, x, q, quantum)
              #  Liqss_reComputeNextTime(Val(1), j, simt, nextStateTime, x, q, quantum)
   
                for k in SD(j)  #j influences k
                    if #= k!=index && =# k!=j  #j will be updated in the next for loop
                        elapsedx = simt - tx[k]; x[k].coeffs[1] = x[k](elapsedx);tx[k] = simt
                        elapsedq = simt - tq[k]; if elapsedq > 0  ;integrateState(Val(1-1),q[k],elapsedq); tq[k] = simt end
                       #=  for b in (jac(k)  )   
                          elapsedq = simt - tq[b]
                          if elapsedq>0 integrateState(Val(1-1),q[b],elapsedq);tq[b]=simt end
                        end   =#
                        clearCache(taylorOpsCache,Val(CS),Val(1));f(k,q,t,taylorOpsCache);computeDerivative(Val(1), x[k], taylorOpsCache[1])
                        Liqss_reComputeNextTime(Val(1), k, simt, nextStateTime, x, q, quantum)
                    end#end if k!=j
                end#end for k depend on j          
              end#end ifcycle check
          # end #end if allow one simulupdate
          end#end if j!=0
        end#end FOR_cycle check
            
      #-------------------------------------------------------------------------------------
      #---------------------------------normal liqss: proceed--------------------------------
      #-------------------------------------------------------------------------------------

        for c in SD(index)   #index influences c       
          elapsedx = simt - tx[c] ;if elapsedx>0 x[c].coeffs[1] = x[c](elapsedx);tx[c] = simt end # 
          elapsedq = simt - tq[c];if elapsedq > 0 ;integrateState(Val(1-1),q[c],elapsedq);tq[c] = simt    end   # c never been visited 
          #= for b in (jac(c)  )    # update other influences
              elapsedq = simt - tq[b] ;if elapsedq>0 integrateState(Val(1-1),q[b],elapsedq);tq[b]=simt  end
          end  =#
          clearCache(taylorOpsCache,Val(CS),Val(1)); f(c,q,t,taylorOpsCache);computeDerivative(Val(1), x[c], taylorOpsCache[1])
          Liqss_reComputeNextTime(Val(1), c, simt, nextStateTime, x, q, quantum)
        end#end for SD
          
    
   
    end#end state/input/event
    #push!(savedVars[index],(x[index][0]+q[index][0])/2)
    push!(savedVars[index],x[index][0])
    push!(savedTimes[index],simt)
 #   push!(savedVarsQ[index],q[index][0])
  
    
  end#end while

#@show temporaryhelper
 #= @timeit "createSol" =# 
 createSol(Val(T),Val(1),savedTimes,savedVars, "singleUpdate$(O)_$(SN)_detection$(Detection)",string(odep.prname),absQ,totalSteps,simulStepCount,0,numSteps,ft)
 #createSol(Val(T),Val(1),savedTimes,savedVars,savedVarsQ, "mliqss$1",string(odep.prname),absQ,totalSteps#= ,totalStepswhenCycles =#,simulStepCount,numSteps,ft)
     # change this to function /constrcutor...remember it is bad to access structs (objects) directly
  
end
