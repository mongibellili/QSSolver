 #using TimerOutputs
function LiQSS_integrate(::Val{O}, s::LiQSS_data{T,Z,O}, odep::NLODEProblem{T,D,Z,Y},f::Function) where {O,T,D,Z,Y}
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
    numSteps = zeros(MVector{T,Int})
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
   # @show a# sys2: [[0.0, 0.01], [-100.0, -100.0]]
    for i = 1:T
      p=1
     
      for k=1:O# deleting this causes scheduler error
        p=p*k
        m=p/k
          u[i][i][k]=p*x[i][k]-m*q[i][k-1]*a[i][i] #  later we will investigate inconsistencies of using data stored vs */ factorial!!! ...also do not confuse getindex for taylor...[0] first element and u[i][1]...first element    ########||||||||||||||||||||||||||||||||||||liqss|||||||||||||||||||||||||||||||||||||||||
          #= println("$p*x[$i][$k]= ", p*x[i][k])
          println("$m*q[$i][$(k-1)]*a[$i][$i]= ",m*q[i][k-1]*a[i][i]) =#
      end
    end
   # @show u #sys2:[[0.0, 19.2], [0.0, 0.0]], [[0.0, 0.0], [1920.0, 0.0]]
    for i = 1:T    
      #= @show u[i][i][2]
      @show x[i][2] =#
      savedVars[i][1].coeffs .= x[i].coeffs  
      quantum[i] = relQ * abs(x[i].coeffs[1]) 
      if quantum[i] < absQ
        quantum[i] = absQ
      end
     # computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
     updateQ(Val(O),i,x,q,quantum,a,u,qaux,olddx,tq,tu,initTime,ft,nextStateTime) 
   #  updateQ(Val(O),i,x,q,quantum,a,u,qaux,olddx,tq,tu,initTime,ft) 
     # computeNextTime(Val(O), i, initTime, nextStateTime, x, quantum)
     #= x[i][1]=a[i][i]*q[i][0]+u[i][i][1]
     x[i][2]=(a[i][i]*q[i][1]+u[i][i][2])/2 =#

    #  Liqss_ComputeNextTime1(Val(O), i, initTime, nextStateTime, x,q, quantum)

     
      clearCache(taylorOpsCache,cacheSize)
      f(i,q,d,t,taylorOpsCache) #+t alloc   change to addT

     
      #@timeit "comp der" 
      computeDerivative(Val(O), x[i], taylorOpsCache[1],integratorCache,0.0)#0.0 used to be elapsed...even down below not neeeded anymore
      #computeDerivative(Val(O), x[j], taylorOpsCache[1])
    #  @timeit "-recompute"
       Liqss_reComputeNextTime(Val(O), i, initTime, nextStateTime, x, q, quantum,a)


      computeNextInputTime(Val(O), i, initTime, 0.1,taylorOpsCache[1] , nextInputTime, x,  quantum)        #   *******************************testing: remove comments later************************
    end
    
    #= @show q[1][0]-x[1][0],q[1][1]-x[1][1],-x[1][2]
    @show q[2][0]-x[2][0],q[2][1]-x[2][1],-x[2][2] =#
    #= @show quantum
    @show x
    @show q
    @show nextStateTime =#
    println("test")

    for i=1:Z
      clearCache(taylorOpsCache,cacheSize)
      output=zcf[i](x,d,t,taylorOpsCache).coeffs[1] #test this evaluation
      oldsignValue[i,2]=output #value
      oldsignValue[i,1]=sign(output) #sign modify 
      computeNextEventTime(i,output,oldsignValue,initTime,  nextEventTime, quantum)#,printCounter)
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
     tempquan=[]
    tempquanTime=[]
    while simt < ft && printcount < 5000000
      printcount+=1
      sch = updateScheduler(nextStateTime,nextEventTime, nextInputTime)
      simt = sch[2]
      if simt==Inf
        println("Equilibrium detected. Simulation ends!")
        count += 1
        
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
         savedTimes[count]=ft
         break
      end
      index = sch[1]
      numSteps[index]+=1
      t[0]=simt
      ##########################################state########################################
      if sch[3] == :ST_STATE
        elapsed = simt - tx[index]
        integrateState(Val(O),x[index],integratorCache,elapsed)#

        tx[index] = simt
        quantum[index] = relQ * abs(x[index].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
        if quantum[index] < absQ
          quantum[index] = absQ
         #=  println("quan limited")
          @show simt =#
        end
       # updateQ(Val(O),index,x,q,quantum,a,u,qaux,olddx,tq,tu,simt,ft)
       
        updateQ(Val(O),index,x,q,quantum,a,u,qaux,olddx,tq,tu,simt,ft,nextStateTime) ########||||||||||||||||||||||||||||||||||||liqss|||||||||||||||||||||||||||||||||||||||||
     
        #= if a[index][index]==0.0
              for b = 1:T # elapsed update all other vars that this derj depends upon.needed for when sys has 3 or more vars.
                #sj = jac[j][b] 
                if jac[index][b]  != 0  
                                                    #   if debug  @show sj   end      
                  elapsedq = simt - tq[b]
                  if elapsedq>0
                  # q[b].coeffs[1] = q[b](elapsed) ## 
                    integrateState(Val(O-1),q[b],integratorCache,elapsedq)
                    tq[b]=simt
                                                    #   if debug println("q$b elapse updated under sj") end
                  end
                end
              end
            
              clearCache(taylorOpsCache,cacheSize)
              f(index,q,d,t,taylorOpsCache)
              #@timeit "comp der" 
              computeDerivative(Val(O), x[index], taylorOpsCache[1],integratorCache,elapsed)
              x2=(2*x[index][2]+olddx[index][2])/2;x1=(x[index][1]+olddx[index][2])/2;xt=x[index][0]
              h=0.0
              if x2!=0.0
              
                      h=sqrt(abs(2*quantum[index]/x2))   #sqrt necessary with u2
                    # q=x+h*x1+h*h*x2/2
                    q[index][0]=xt-h*h*x2/2
                    q[index][1]=x1+h*x2
                    # q=x-h*h*u2/2
                    # q1=u1+h*u2  #(250 allocations: 16.07 KiB)
                    #  @show h,x2,i,u1,u2
              else
                      println("x2==0")
                      if x1!=0.0
                          h=abs(quantum[index]/x1)
                          q[index][0]=xt+h*x1
                          q[index][1]=x1
                      else
                        q[index][0]=xt
                        q[index][1]=x1
                      end
                  #=q1=x1+h*x2  #(252 allocations: 16.07 KiB)
                  else
                      q1=x1 =#
                      #println("ddx=0")
              end 
              if abs(q[index][0]-xt)>quantum[index]+1e-12#uncommenting this did nothing
                  println("at end of q update: else side, q-x >quan !!!! simt= ",simt,", q-x= ",q[index][0]-xt,", quan= ",quantum[index])
              end
              nextStateTime[index]=simt+h

        end =#
        tq[index]=simt
        # @timeit "computeNextTime" 
     
       # computeNextTime(Val(O), index, simt, nextStateTime, x, quantum) #
       # Liqss_ComputeNextTime(Val(O), index, simt, nextStateTime, x,q, quantum)
        for i = 1:length(SD[index])
          j = SD[index][i] 
          if j != 0           
            elapsedx = simt - tx[j]
            if elapsed > 0
              #"evaluate" x at new time only...derivatives get updated next using computeDerivativ()
              x[j].coeffs[1] = x[j](elapsedx)
            #=  quantum[j] = relQ * abs(x[j].coeffs[1])  #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
              if quantum[j] < absQ
                quantum[j] = absQ
                println("quan limited inside dependency")
                @show simt
              end =#
             # q[j].coeffs[1] = q[j](elapsed)
              tx[j] = simt
             # tq[j] = simt
            end
            elapsedq = simt - tq[j]
            if elapsedq > 0
              #"evaluate" x only at new time ...derivatives get updated next using computeDerivativ()
              
              integrateState(Val(O-1),q[j],integratorCache,elapsedq)
             # q[j].coeffs[1] = q[j](elapsedq) # ouch ! this bit me for a day: q needs to be updated here for recomputeNext, not just for the derivatives!!!
                                                      # if debug println("x$j elapse updated")  end 
              
             tq[j] = simt
            end
            

            for b = 1:T # elapsed update all other vars that this derj depends upon.needed for when sys has 3 or more vars.
              #sj = jac[j][b] 
              if jac[j][b]  != 0  
                                                  #   if debug  @show sj   end      
                elapsedq = simt - tq[b]
                if elapsedq>0
                 # q[b].coeffs[1] = q[b](elapsed) ## 
                  integrateState(Val(O-1),q[b],integratorCache,elapsedq)
                  tq[b]=simt
                                                   #   if debug println("q$b elapse updated under sj") end
                end
              end
            end
           
            clearCache(taylorOpsCache,cacheSize)
            f(j,q,d,t,taylorOpsCache)
            #@timeit "comp der" 
          #  if x[j][1]!=taylorOpsCache[1][0]
            computeDerivative(Val(O), x[j], taylorOpsCache[1],integratorCache,elapsed)
            #computeDerivative(Val(O), x[j], taylorOpsCache[1])
          #  @timeit "-recompute"
             Liqss_reComputeNextTime(Val(O), j, simt, nextStateTime, x, q, quantum,a)
          #  end
             #= if 50>simt>9
              #@show a[1][2], a[2][1]
             # @show u[j][index][1]+a[j][j]*q[j][0]+a[j][index]*q[index][0]
             # @show x[j][1]
             # @show nextStateTime
              @show simt
              @show printcount
              #limitedPrint-=1
             end =#
          end#end if j!=0
        end#end for SD
        for i = 1:length(SZ[index])
          j = SZ[index][i] 
          if j != 0             
            #normally and later i should update q (integrate q=q+e derQ  for higher orders)
            clearCache(taylorOpsCache,cacheSize)
            computeNextEventTime(j,zcf[j](x,d,t,taylorOpsCache)[0],oldsignValue,simt,  nextEventTime, quantum)#,printCounter)
          end  #end if j!=0
        end#end for SZ

       # @timeit "updateLinearApprox" 
        updateLinearApprox(Val(O),index,x,q,a,u,qaux,olddx,tu,simt)########||||||||||||||||||||||||||||||||||||liqss|||||||||||||||||||||||||||||||||||||||||
        ##################################input########################################
      elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself  
       # println("input")
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
       # println("ebent")
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

   #=    @show x
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
         if index==2
          push!(tempquanTime, simt)
          push!(tempquan, x[index][0]-q[index][0])
        end
      end#end if save
    end#end while
    for i=1:T# throw away empty points
      resize!(savedVars[i],count)
    end
  #  print_timer()
   @show printcount
   println("------------------------------------------------------------------------------------")
   #=  display(plot!(tempquanTime, tempquan,title="quantum graph for x2",label="liqss$O")) 
  println("press enter to exit")
  readline() =#
    resize!(savedTimes,count)
    Sol(savedTimes, savedVars,"liqss$O",string(nameof(f)),numSteps,absQ)
    end#end integrate
    
    