#using TimerOutputs
function QSS_integrate(CommonqssData::CommonQSS_data{O,Z}, odep::NLODEProblem{PRTYPE,T,Z,Y,CS},f::Function,jac::Function,SD::Function,map::Function) where {PRTYPE,O,T,Z,Y,CS}

  #cacheA=specialLiqssData.cacheA
  ft = CommonqssData.finalTime;initTime = CommonqssData.initialTime;relQ = CommonqssData.dQrel;absQ = CommonqssData.dQmin;maxErr=CommonqssData.maxErr;

  #savetimeincrement=CommonqssData.savetimeincrement;savetime = savetimeincrement
  quantum = CommonqssData.quantum;nextStateTime = CommonqssData.nextStateTime;nextEventTime = CommonqssData.nextEventTime;nextInputTime = CommonqssData.nextInputTime
  tx = CommonqssData.tx;tq = CommonqssData.tq;x = CommonqssData.x;q = CommonqssData.q;t=CommonqssData.t
  savedVars=CommonqssData.savedVars;
  savedTimes=CommonqssData.savedTimes;integratorCache=CommonqssData.integratorCache;taylorOpsCache=CommonqssData.taylorOpsCache;#cacheSize=odep.cacheSize
  #prevStepVal = specialLiqssData.prevStepVal
  #*********************************problem info*****************************************
  d = odep.discreteVars
  

  zc_SimpleJac=odep.ZCjac

  HZ=odep.HZ
  HD=odep.HD
  SZ=odep.SZ
 
  evDep = odep.eventDependencies

  if VERBOSE @show HD,HZ end
 

  #qaux=liqssdata.qaux;dxaux=liqssdata.dxaux#= olddx=liqssdata.olddx; ; olddxSpec=liqssdata.olddxSpec =#


 #=  pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
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
  exacteA(q,cacheA,1,1) =#
  #@show exacteA
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
  computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
  #updateQ(Val(O),i,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,initTime,ft,nextStateTime) 
end
#= for i = 1:T
   clearCache(taylorOpsCache,Val(CS),Val(O));f(i,-1,-1,q,d,t,taylorOpsCache);
   computeDerivative(Val(O), x[i], taylorOpsCache[1])#0.0 used to be elapsed...even down below not neeeded anymore
  Liqss_reComputeNextTime(Val(O), i, initTime, nextStateTime, x, q, quantum)
  computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)#not complete, currently elapsed=0.1 is temp until fixed
 #prevStepVal[i]=x[i][0]#assignXPrevStepVals(Val(O),prevStepVal,x,i)
end =#
  

  #assignXPrevStepVals(Val(O),prevStepVal,x,i)
  


#@show nextStateTime,nextEventTime
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
simt = initTime ;totalSteps=0;prevStepTime=initTime;modifiedIndex=0; countEvents=0
  
while simt<ft && totalSteps < 50000000
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
    numSteps[index]+=1;
    elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
    quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
    for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt    
    computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) #
    for j in (SD(index))
        elapsedx = simt - tx[j];if elapsedx > 0 x[j].coeffs[1] = x[j](elapsedx);tx[j] = simt end
        # quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]         
        elapsedq = simt - tq[j];if elapsedq > 0 integrateState(Val(O-1),q[j],elapsedq);tq[j] = simt  end#q needs to be updated here for recomputeNext        
        for b in (jac(j)  )    
          elapsedq = simt - tq[b]
          if elapsedq>0
            integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt
          end
         end
        clearCache(taylorOpsCache,Val(CS),Val(O));f(j,-1,-1,q,d,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])  
        reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
    end#end for SD

    for j in (SZ[index])
      for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
          #elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
      end            
      clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,q,d,t,taylorOpsCache)   # run ZCF--------      
      computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
  end#end for SZ
  #=  if 8.377869511741662<=simt<=14.578335986845602 && (index==2 #= || index==5 =#) 
      println("intgrator at simt=$simt index=$index")
     
      @show x[index],nextStateTime[index]
      
     end =#

    ##################################input########################################
  elseif stepType == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
  #=   if VERBOSE println("nmliqss discreteintgrator under input, index= $index, totalsteps= $totalSteps")  end
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
       # quantum[j] = relQ * abs(x[j].coeffs[1]) ;quantum[j]=quantum[j] < absQ ? absQ : quantum[j];quantum[j]=quantum[j] > maxErr ? maxErr : quantum[j]   
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
             
            elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
           # elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
       
        end              
       #=  clearCache(taylorOpsCache,Val(CS),Val(O))#normally and later i should update x,q (integrate q=q+e derQ  for higher orders)
        computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer)  =#
        clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,x,d,t,taylorOpsCache)        
         computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
     
    end =#
  #################################################################event########################################
  else
         # printcounter=5
         if VERBOSE println("x at start of event simt=$simt index=$index") end
       
       
      
          for b in zc_SimpleJac[index] # elapsed update all other vars that this zc depends upon.
              
             #=  elapsedq = simt - tq[b];
              if elapsedq>0 
                integrateState(Val(O),q[b],elapsedq);
               # println("should update var1 b=",b)
                tq[b]=simt 
              end =#
             elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
           
          end    
         # modifiedIndex=0#first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in O to know ev+ or ev- 
        
          clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,index,-1,q,d,t,taylorOpsCache)    # run ZCF-------- 
          if VERBOSE  @show oldsignValue[index,2],taylorOpsCache[1][0]  end
          #=  println(" just after event")
        
          @show t 
          @show nextEventTime   =# 
                
         #=  if (taylorOpsCache[1][0])>-1e-13  #zcf>=0      # sign is not needed here
            modifiedIndex=2*index-1   # the  event that just occured is at  this index
            #println("event1")
          else
            modifiedIndex=2*index
          end  =# 
          if oldsignValue[index,2]*taylorOpsCache[1][0]>0 && abs(taylorOpsCache[1][0])>1e-13 # if both have same sign and zcf is not very small
            computeNextEventTime(index,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum) 
            println("case wrong estimation of event at $simt")
            continue
          end

         #=  if oldsignValue[index,2]<0.0 #rise
            modifiedIndex=2*index-1 
          else # drop
            modifiedIndex=2*index

          end =#
          if taylorOpsCache[1][0]>oldsignValue[index,2] #scheduled rise
            modifiedIndex=2*index-1 
          elseif taylorOpsCache[1][0]<oldsignValue[index,2] #scheduled drop
            modifiedIndex=2*index
          else # == coming from a triggered zc (very rare when a check coincides with zcf==0)
            #= if oldsignValue[index,2]<0.0 #rise
              modifiedIndex=2*index-1 
            else
              modifiedIndex=2*index
            end =#
          end
          oldsignValue[index,2]=taylorOpsCache[1][0]
          oldsignValue[index,1]=sign(taylorOpsCache[1][0])
        # println("x before event= ",x) 
        if VERBOSE @show modifiedIndex end
           
         for b in evDep[modifiedIndex].evContRHS
                
                  elapsedq = simt - tq[b];
                  if elapsedq>0 
                    integrateState(Val(O-1),q[b],elapsedq);
                  #  println("var2 here b = ",b)
                    
                    tq[b]=simt 
                   end
               
           end
           
         ###### eventf[modifiedIndex](x,d,t,taylorOpsCache) #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
        # @show x,modifiedIndex
         f(-1,-1,modifiedIndex,x,d,t,taylorOpsCache)# execute event--------
      
          for i in evDep[modifiedIndex].evCont
            #------------event influences a Continete var: x already updated in event...here update quantum and q and computenext
            
                quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
                q[i][0]=x[i][0];tx[i] = simt;tq[i] = simt # for liqss updateQ?
             #   firstguess=updateQ(Val(O),i,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,simt,ft,nextStateTime) ;tq[i] = simt   
                computeNextTime(Val(O), i, simt, nextStateTime, x, quantum) 
           
          end
          nextEventTime[index]=Inf   #investigate more 
         # computeNextEventTime(index,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum) # it could happen zcf=0.0 then infinite event

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
              @show simt,j,x
           
          end
          for j in (HZ[modifiedIndex])
                
                    for b in zc_SimpleJac[j] # elapsed update all other vars that this derj depends upon.
                          
                        #elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],elapsedx);tx[b]=simt end
                       elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
                    
                    end            
                  #= clearCache(taylorOpsCache,Val(CS),Val(O)) #normally and later i should update q (integrate q=q+e derQ  for higher orders)          
                  computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer) =#
                  
                  clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,q,d,t,taylorOpsCache)  # run ZCF-------- 
                  if VERBOSE @show j,oldsignValue[j,2],taylorOpsCache[1][0] end     
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
    
         #=  for i =1:T 
            push!(savedVars[i],x[i][0])
            push!(savedTimes[i],simt)
          end =#
      else
        countEvents+=1
        # if 1e-3<simt<1e-2 || 1e-5<simt<1e-4
        for j in (HD[modifiedIndex])
          
          push!(savedVars[j],x[j][0])
          push!(savedTimes[j],simt)
          
        end
      # end
      end
      
    #push!(savedVarsQ[i],q[i][0])
    prevStepTime=simt
# end
end#end while
 
@show countEvents
#@show savedVars
#createSol(Val(T),Val(O),savedTimes,savedVars, "qss$O",string(nameof(f)),absQ,totalSteps,0)#0 I track simulSteps 
createSol(Val(T),Val(O),savedTimes,savedVars, "qss$O",string(odep.prname),absQ,totalSteps,0,numSteps,ft)
end#end integrate