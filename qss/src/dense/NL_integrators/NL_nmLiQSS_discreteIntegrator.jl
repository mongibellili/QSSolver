#using TimerOutputs
function nmLiQSS_integrate(CommonqssData::CommonQSS_data{O,Z},liqssdata::LiQSS_data{O,false},specialLiqssData::SpecialLiqssQSS_data, odep::NLODEProblem{PRTYPE,T,Z,Y,CS},f::Function,jac::Function,SD::Function,exacteA::Function) where {PRTYPE,O,T,Z,Y,CS}
  cacheA=specialLiqssData.cacheA
  ft = CommonqssData.finalTime;initTime = CommonqssData.initialTime;relQ = CommonqssData.dQrel;absQ = CommonqssData.dQmin;maxErr=CommonqssData.maxErr;
  savetimeincrement=CommonqssData.savetimeincrement;savetime = savetimeincrement
  quantum = CommonqssData.quantum;nextStateTime = CommonqssData.nextStateTime;nextEventTime = CommonqssData.nextEventTime;nextInputTime = CommonqssData.nextInputTime
  tx = CommonqssData.tx;tq = CommonqssData.tq;x = CommonqssData.x;q = CommonqssData.q;t=CommonqssData.t
  savedVars=CommonqssData.savedVars;
  savedTimes=CommonqssData.savedTimes;integratorCache=CommonqssData.integratorCache;taylorOpsCache=CommonqssData.taylorOpsCache;cacheSize=odep.cacheSize
  #prevStepVal = specialLiqssData.prevStepVal
  #*********************************problem info*****************************************
  d = odep.discreteVars
  

  zc_SimpleJac=odep.ZCjac

  HZ=odep.HZ
  HD=odep.HD
  SZ=odep.SZ
 
  evDep = odep.eventDependencies



  qaux=liqssdata.qaux;dxaux=liqssdata.dxaux#= olddx=liqssdata.olddx; ; olddxSpec=liqssdata.olddxSpec =#


  pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
 # savedVarsQ = Vector{Vector{Float64}}(undef, T)  

                            
  respp = pointer(Vector{Float64}(undef, 2))
  trackSimul = Vector{Int}(undef, 1)
 # cacheRatio=zeros(5);cacheQ=zeros(5)
  cacherealPosi=Vector{Vector{Float64}}(undef,3);cacherealPosj=Vector{Vector{Float64}}(undef,3);
 #=  for i=1:T
    savedVarsQ[i]=Vector{Float64}()     
  end =#
for i =1:3
  cacherealPosi[i]=zeros(2)
  cacherealPosj[i]=zeros(2)
end
  exacteA(q,cacheA,1,1)

 # @show HD,HZ,SZ,d
 #@show f
  #********************************helper values*******************************  

  oldsignValue = MMatrix{Z,2}(zeros(Z*2))  #usedto track if zc changed sign; each zc has a value and a sign 
  numSteps = Vector{Int}(undef, T)
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
  numSteps[i]=0
   push!(savedVars[i],x[i][0])
   push!(savedTimes[i],0.0)
  
  quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
  
  updateQ(Val(O),i,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,initTime,ft,nextStateTime) 
end
for i = 1:T
   clearCache(taylorOpsCache,Val(CS),Val(O));f(i,-1,-1,q,d,t,taylorOpsCache);
   computeDerivative(Val(O), x[i], taylorOpsCache[1])#0.0 used to be elapsed...even down below not neeeded anymore
  Liqss_reComputeNextTime(Val(O), i, initTime, nextStateTime, x, q, quantum)
  computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)#not complete, currently elapsed=0.1 is temp until fixed
 #prevStepVal[i]=x[i][0]#assignXPrevStepVals(Val(O),prevStepVal,x,i)
end
  

  #assignXPrevStepVals(Val(O),prevStepVal,x,i)
  


#@show nextStateTime,nextInputTime
for i=1:Z
  #= clearCache(taylorOpsCache,Val(CS),Val(O));output=zcf[i](x,d,t,taylorOpsCache).coeffs[1]  =#
  clearCache(taylorOpsCache,Val(CS),Val(O));
  #@show taylorOpsCache
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
simt = initTime ;totalSteps=0;prevStepTime=initTime;modifiedIndex=0
  
while simt < ft && totalSteps < 50000000
  sch = updateScheduler(Val(T),nextStateTime,nextEventTime, nextInputTime)
  simt = sch[2];index = sch[1];stepType=sch[3]
 # @timeit "saveLast" 
   if  simt>ft  
    #saveLast!(Val(T),Val(O),savedVars, savedTimes,saveVarsHelper,ft,prevStepTime, x)
    break   ###################################################break##########################################
  end
  totalSteps+=1

  t[0]=simt
  ##########################################state######################################## 
  if stepType == :ST_STATE
    xitemp=x[index][0]
    numSteps[index]+=1;
    elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    dirI=x[index][0]-xitemp
    for b in (jac(index)  )    # update Qb : to be used to calculate exacte Aindexb
      elapsedq = simt - tq[b] ;
      if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
    end
    firstguess=updateQ(Val(O),index,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,simt,ft,nextStateTime) ;tq[index] = simt   
    #----------------------------------------------------check dependecy cycles---------------------------------------------   
    trackSimul[1]=0 
    #= for i =1:5
      cacheRatio[i]=0.0; cacheQ[i]=0.0; 
    end       =#       
    for j in SD(index)
      for b in (jac(j)  )    # update Qb: to be used to calculate exacte Ajb
        elapsedq = simt - tq[b] ;
        if elapsedq>0  integrateState(Val(O-1),q[b],elapsedq); tq[b]=simt  end
      end
      exacteA(q,cacheA,index,j);aij=cacheA[1]# can be passed to simul so that i dont call exactfunc again
      exacteA(q,cacheA,j,index);aji=cacheA[1]
     
    
    
      if j!=index && aij*aji!=0.0
          #prvStepValj= savedVars[j][end]#getPrevStepVal(prevStepVal,j) 
          for i =1:3
            cacherealPosi[i][1]=0.0; cacherealPosi[i][2]=0.0
            cacherealPosj[i][1]=0.0; cacherealPosj[i][2]=0.0
          end 
          if nmisCycle_and_simulUpdate(cacherealPosi,cacherealPosj,aij,aji,respp,pp,trackSimul,Val(O),index,j,dirI,firstguess,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,simt,ft)
            simulStepCount+=1
           clearCache(taylorOpsCache,Val(CS),Val(O));f(index,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[index], taylorOpsCache[1])
          #  clearCache(taylorOpsCache,Val(CS),Val(O));f(j,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
          # Liqss_reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum)
          #  Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)

            for k in SD(j)  #j influences k
                if k!=index && k!=j
                    elapsedx = simt - tx[k]; x[k].coeffs[1] = x[k](elapsedx);tx[k] = simt
                    elapsedq = simt - tq[k]; if elapsedq > 0  ;integrateState(Val(O-1),q[k],elapsedq); tq[k] = simt end
                    for b in (jac(k)  )   
                      elapsedq = simt - tq[b]
                      if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
                    end  
                    clearCache(taylorOpsCache,Val(CS),Val(O));f(k-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[k], taylorOpsCache[1])
                    Liqss_reComputeNextTime(Val(O), k, simt, nextStateTime, x, q, quantum)
                end#end if k!=0
            end#end for k depend on j     
            
            
          end#end ifcycle check
      # end #end if allow one simulupdate
      end#end if j!=0
    end#end FOR_cycle check
        

  if trackSimul[1]!=0  #qi changed after throw
   # clearCache(taylorOpsCache,Val(CS),Val(O));f(index,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[index], taylorOpsCache[1])
    Liqss_reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum)
  end

  #-------------------------------------------------------------------------------------
  #---------------------------------normal liqss: proceed--------------------------------
  #-------------------------------------------------------------------------------------

    for c in SD(index)   #index influences c       
      elapsedx = simt - tx[c] ;
      if elapsedx>0 
        x[c].coeffs[1] = x[c](elapsedx);
        tx[c] = simt

       end # 
      elapsedq = simt - tq[c];if elapsedq > 0 ;integrateState(Val(O-1),q[c],elapsedq);tq[c] = simt    end   # c never been visited 
      #= for b in (jac(c)  )    # update other influences
          elapsedq = simt - tq[b] ;if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt  end
      end  =#
      clearCache(taylorOpsCache,Val(CS),Val(O)); f(c,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[c], taylorOpsCache[1])
      Liqss_reComputeNextTime(Val(O), c, simt, nextStateTime, x, q, quantum)
    end#end for SD

    for j in (SZ[index])
      for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],elapsedx);tx[b]=simt end
          #elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
      end            
      clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,x,d,t,taylorOpsCache)        
      computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
  end#end for SZ
  #=  if 8.377869511741662<=simt<=14.578335986845602 && (index==2 #= || index==5 =#) 
      println("intgrator at simt=$simt index=$index")
     
      @show x[index],nextStateTime[index]
      
     end =#

    ##################################input########################################
  elseif stepType == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
    
    elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt 
      for b in jac(index) 
        elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
      end
    clearCache(taylorOpsCache,Val(CS),Val(O));f(index,-1,-1,q,d,t,taylorOpsCache)
    computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
    computeDerivative(Val(O), x[index], taylorOpsCache[1])

    for j in(SD(index))  
      elapsedx = simt - tx[j];
      if elapsedx > 0 
        x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt 
        quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]   
      end
      elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext                 
      # elapsed update all other vars that this derj depends upon.
        for b in jac(j) 
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
        end
      
        clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
    end#end for
    for j in (SZ[index])
      
        for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
             
            elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],elapsedx);tx[b]=simt end
           # elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
       
        end              
       #=  clearCache(taylorOpsCache,Val(CS),Val(O))#normally and later i should update x,q (integrate q=q+e derQ  for higher orders)
        computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer)  =#
        clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,x,d,t,taylorOpsCache)        
         computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
     
    end
  #################################################################event########################################
  else
         # printcounter=5
         #=  println("x at start of event")
          @show x
          @show simt 
          @show quantum =#
      
          for b in zc_SimpleJac[index] # elapsed update all other vars that this zc depends upon.
              
              elapsedx = simt - tx[b];
              if elapsedx>0 
                integrateState(Val(O),x[b],elapsedx);
               # println("should update var1 b=",b)
                tx[b]=simt 
              end
             elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
           
          end    
         # modifiedIndex=0#first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 
         
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
                    integrateState(Val(O),x[b],elapsedx);
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
              elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt;#= @show q[j] =#  end#q needs to be updated here for recomputeNext                 
              for b = 1:T # elapsed update all other vars that this derj depends upon.
                if b in jac(j)   
                  elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt;#= @show q[b] =# end
                end
              end
              clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
              reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
              #@show j,x
           
          end
          for j in (HZ[modifiedIndex])
                 
                    for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
                          
                        elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],elapsedx);tx[b]=simt end
                       #elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
                    
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
  #for i=1:T
  if stepType != :ST_EVENT
    push!(savedVars[index],x[index][0])
    push!(savedTimes[index],simt)
  else
    for j in (HD[modifiedIndex])
      push!(savedVars[j],x[j][0])
      push!(savedTimes[j],simt)
    end
  end
    #push!(savedVarsQ[i],q[i][0])
# end
end#end while
 

#createSol(Val(T),Val(O),savedTimes,savedVars, "qss$O",string(nameof(f)),absQ,totalSteps,0)#0 I track simulSteps 
createSol(Val(T),Val(O),savedTimes,savedVars, "nmLiqss$O",string(odep.prname),absQ,totalSteps,0,numSteps,ft)
end#end integrate