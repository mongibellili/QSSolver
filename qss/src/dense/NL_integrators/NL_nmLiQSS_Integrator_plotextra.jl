
function nmLiQSS_integrate(CommonqssData::CommonQSS_data{O,0},liqssdata::LiQSS_data{O,false},specialLiqssData::SpecialLiqssQSS_data, odep::NLODEProblem{PRTYPE,T,0,0,CS},f::Function,jac::Function,SD::Function,exacteA::Function) where {PRTYPE,CS,O,T}
  cacheA=specialLiqssData.cacheA
  direction=specialLiqssData.direction
  #qminus= specialLiqssData.qminus
  buddySimul=specialLiqssData.buddySimul
  ft = CommonqssData.finalTime;initTime = CommonqssData.initialTime;relQ = CommonqssData.dQrel;absQ = CommonqssData.dQmin;maxErr=CommonqssData.maxErr;
  
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
  simulStepsVals = Vector{Vector{Float64}}(undef, T)
  simulStepsDers = Vector{Vector{Float64}}(undef, T)
  simulStepsTimes = Vector{Vector{Float64}}(undef, T)
  simulHTimes=Vector{Float64}()
    simulHVals=Vector{Float64}()
    simulqxiVals=Vector{Float64}()
    simulqxjVals=Vector{Float64}()
    simuldeltaiVals=Vector{Float64}()
    simuldeltajVals=Vector{Float64}()
  exacteA(q,cacheA,1,1)  # this 'unnecessary call' 'compiles' the function and it helps remove allocations when used after !!!

  #@show exacteA
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
    numSteps[i]=0
    simulStepsTimes[i]=Vector{Float64}()
    simulStepsVals[i]=Vector{Float64}()
    simulStepsDers[i]=Vector{Float64}()
    #= @timeit "savevars" =# push!(savedVars[i],x[i][0])
     push!(savedTimes[i],0.0)
     quantum[i] = relQ * abs(x[i].coeffs[1]) ;quantum[i]=quantum[i] < absQ ? absQ : quantum[i];quantum[i]=quantum[i] > maxErr ? maxErr : quantum[i] 
    updateQ(Val(O),i,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,initTime,ft,nextStateTime) 
  end
  for i = 1:T
     clearCache(taylorOpsCache,Val(CS),Val(O));f(i,q,t,taylorOpsCache);
     computeDerivative(Val(O), x[i], taylorOpsCache[1])#0.0 used to be elapsed...even down below not neeeded anymore
    Liqss_reComputeNextTime(Val(O), i, initTime, nextStateTime, x, q, quantum)
    computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)#not complete, currently elapsed=0.1 is temp until fixed
   #prevStepVal[i]=x[i][0]#assignXPrevStepVals(Val(O),prevStepVal,x,i)
  end


  ###################################################################################################################################################################
  ####################################################################################################################################################################
  #---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
  ###################################################################################################################################################################
  #################################################################################################################################################################### 
  simt = initTime ;simulStepCount=0;totalSteps=0;
 
  #simul=false
breakcounter=3

  while simt < ft && totalSteps < 300000000
    sch = updateScheduler(Val(T),nextStateTime,nextEventTime, nextInputTime)
    simt = sch[2];index = sch[1]
    if simt>=ft# this is needed so that interpolated solution be defined...can rewrite 'interpolated' and remove this check. also in updateQ: divide by 1-ha 
      #simt=ft
      break
    end
    if breakcounter<0
      break
    end
    numSteps[index]+=1;totalSteps+=1
    t[0]=simt
    ##########################################state########################################
    if sch[3] == :ST_STATE
      xitemp=x[index][0]
        elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
        quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index] 
       # dirI=x[index][0]-savedVars[index][end]  
        dirI=x[index][0]-xitemp
        for b in (jac(index)  )    # update Qb : to be used to calculate exacte Aindexb
          elapsedq = simt - tq[b] ;
          if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
        end
        firstguessH=updateQ(Val(O),index,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,simt,ft,nextStateTime) ;tq[index] = simt   
        #----------------------------------------------------check dependecy cycles---------------------------------------------                
        for j in SD(index)
          for b in (jac(j)  )    # update Qb: to be used to calculate exacte Ajb
            elapsedq = simt - tq[b] ;
            if elapsedq>0  integrateState(Val(O-1),q[b],elapsedq); tq[b]=simt  end
          end
          exacteA(q,cacheA,index,j);aij=cacheA[1]
          exacteA(q,cacheA,j,index);aji=cacheA[1]
          exacteA(q,cacheA,index,index);aii=cacheA[1]
          exacteA(q,cacheA,j,j);ajj=cacheA[1]
         #=  exacteA(x,cacheA,index,j);aij=cacheA[1]
          exacteA(x,cacheA,j,index);aji=cacheA[1] =#
        
        
          if j!=index && aij*aji!=0.0 #= && 1.2*abs(aij*aji)>abs(aii*ajj) =#
              #prvStepValj= #getPrevStepVal(prevStepVal,j)  
              if nmisCycle_and_simulUpdate(simuldeltaiVals,simuldeltajVals, simulqxiVals,simulqxjVals,  simulHTimes,simulHVals,Val(O),index,j,dirI,firstguessH,x,q,quantum,exacteA,cacheA,dxaux,qaux,tx,tq,simt,ft)
                simulStepCount+=1
                push!(simulStepsTimes[index],simt)
                push!(simulStepsVals[index],x[index][0])
                push!(simulStepsDers[index],x[index][1])
                push!(simulStepsTimes[j],simt)
                push!(simulStepsVals[j],x[j][0])
                push!(simulStepsDers[j],x[j][1])
               # breakcounter-=1
                clearCache(taylorOpsCache,Val(CS),Val(O));f(index,q,t,taylorOpsCache);computeDerivative(Val(O), x[index], taylorOpsCache[1])
              #  clearCache(taylorOpsCache,Val(CS),Val(O));f(j,q,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1])
                Liqss_reComputeNextTime(Val(O), index, simt, nextStateTime, x, q, quantum)
              #  Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
   
                for k in SD(j)  #j influences k
                    if k!=index && k!=j
                        elapsedx = simt - tx[k]; x[k].coeffs[1] = x[k](elapsedx);tx[k] = simt
                        elapsedq = simt - tq[k]; if elapsedq > 0  ;integrateState(Val(O-1),q[k],elapsedq); tq[k] = simt end
                        for b in (jac(k)  )   
                          elapsedq = simt - tq[b]
                          if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
                        end  
                        clearCache(taylorOpsCache,Val(CS),Val(O));f(k,q,t,taylorOpsCache);computeDerivative(Val(O), x[k], taylorOpsCache[1])
                        Liqss_reComputeNextTime(Val(O), k, simt, nextStateTime, x, q, quantum)
                    end#end if k!=0
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
          elapsedq = simt - tq[c];if elapsedq > 0 ;integrateState(Val(O-1),q[c],elapsedq);tq[c] = simt    end   # c never been visited 
          #= for b in (jac(c)  )    # update other influences
              elapsedq = simt - tq[b] ;if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt  end
          end  =#
          clearCache(taylorOpsCache,Val(CS),Val(O)); f(c,q,t,taylorOpsCache);computeDerivative(Val(O), x[c], taylorOpsCache[1])
          Liqss_reComputeNextTime(Val(O), c, simt, nextStateTime, x, q, quantum)
        end#end for SD
          
     
      ##################################input########################################
    elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself    
     @show index
      elapsed = simt - tx[index];integrateState(Val(O),x[index],elapsed);tx[index] = simt 
      quantum[index] = relQ * abs(x[index].coeffs[1]) ;quantum[index]=quantum[index] < absQ ? absQ : quantum[index];quantum[index]=quantum[index] > maxErr ? maxErr : quantum[index]   
      for k = 1:O q[index].coeffs[k] = x[index].coeffs[k] end; tq[index] = simt 
        for b in jac(index) 
          elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
        end
       clearCache(taylorOpsCache,Val(CS),Val(O));f(index,q,t,taylorOpsCache)
      computeNextInputTime(Val(O), index, simt, elapsed,taylorOpsCache[1] , nextInputTime, x,  quantum)
      computeDerivative(Val(O), x[index], taylorOpsCache[1]#= ,integratorCache,elapsed =#)
  
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
           clearCache(taylorOpsCache,Val(CS),Val(O));f(j,q,t,taylorOpsCache);computeDerivative(Val(O), x[j], taylorOpsCache[1]#= ,integratorCache,elapsed =#)
          reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum)
      end#end for
     #=  for i = 1:length(SZ[index])
        j = SZ[index][i] 
        if j != 0   
          for b = 1:T # elapsed update all other vars that this derj depends upon.
            if zc_SimpleJac[j][b] != 0     
              elapsedx = simt - tx[b];if elapsedx>0 integrateState(Val(O),x[b],integratorCache,elapsedx);tx[b]=simt end
             # elapsedq = simt - tq[b];if elapsedq>0 integrateState(Val(O-1),q[b],elapsedq);tq[b]=simt end
            end
          end              
         #=   clearCache(taylorOpsCache,Val(CS),Val(O))#normally and later i should update x,q (integrate q=q+e derQ  for higher orders)
          computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)#,maxIterer)  =#
           clearCache(taylorOpsCache,Val(CS),Val(O));f(-1,j,-1,x,d,t,taylorOpsCache)        
                   computeNextEventTime(j,taylorOpsCache[1],oldsignValue,simt,  nextEventTime, quantum)
        end  
      end =#
    #################################################################event########################################
    end#end state/input/event
    push!(savedVars[index],x[index][0])
    push!(savedTimes[index],simt)
  end#end while
 #=  p1=plot(simulHTimes,simulHVals,marker=(:circle),markersize=2,ylims=(0.0,0.2))
  savefig(p1, "_tyson___H01bothanaly_v0205.png") =#
 #=  p1i=plot()
  p1i=plot!(p1i,simulHTimes,simulqxiVals,marker=(:circle),markersize=2,ylims=(0.0,0.01),label="|x-q|")
  p1i=plot!(p1i,simulHTimes,simuldeltaiVals,marker=(:star),markersize=2,linestyle=:dash,ylims=(0.0,0.01),label="delta")
  savefig(p1i, "_tyson___qxiitertotalv0201.png")
  p1j=plot()
  p1j=plot!(p1j,simulHTimes,simulqxjVals,marker=(:circle),markersize=2,ylims=(0.0,0.01),label="|x-q|")
  p1j=plot!(p1j,simulHTimes,simuldeltajVals,marker=(:star),markersize=2,linestyle=:dash,ylims=(0.0,0.01),label="delta")
  savefig(p1j, "_tyson___qxjitertotalv0201.png") =#

 #=  p1i=plot()
  p1i=plot!(p1i,simulHTimes,simulqxiVals,marker=(:circle),markersize=2,ylims=(0.0,2*1e-5),label="|x-q|")
  p1i=plot!(p1i,simulHTimes,simuldeltaiVals,marker=(:star),markersize=2,linestyle=:dash,ylims=(0.0,2*1e-5),label="delta")
  savefig(p1i, "_tyson___qxiibothanaly_v021e_6.png")
  p1j=plot()
  p1j=plot!(p1j,simulHTimes,simulqxjVals,marker=(:circle),markersize=2,ylims=(0.0,2*1e-5),label="|x-q|")
  p1j=plot!(p1j,simulHTimes,simuldeltajVals,marker=(:star),markersize=2,linestyle=:dash,ylims=(0.0,2*1e-5),label="delta")
  savefig(p1j, "_tyson___qxjibothanaly_v021e_6.png") =#
 #=  readline()
  println("press keyboard") =#

 #= @timeit "createSol" =# 
 createSol(Val(T),Val(O),savedTimes,savedVars, "nmliqss$O",string(odep.prname),absQ,totalSteps,simulStepCount,numSteps,ft,simulStepsVals,simulStepsDers,simulStepsTimes)
     # change this to function /constrcutor...remember it is bad to access structs (objects) directly
  
end
