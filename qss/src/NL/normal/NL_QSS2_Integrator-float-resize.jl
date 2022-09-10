function QSS_integrate(::Val{2}, s::QSS_data, settings::SimSetting, odep::NLODEProblem)
 
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
x = s.x
q = s.q
#=    x=Taylor0{Float64}[]
q=Taylor0{Float64}[] =#
#*********************************problem info*****************************************
initConditions = odep.initConditions
d = odep.discreteVars
#----------to compute derivatives
jacobian = odep.jacobian
discJac = odep.discreteJacobian
#----------to compute ZC expressions
zc_jac = odep.ZC_jacobian
ZC_discJac = odep.ZC_jacDiscrete
#-----------to execute event Dependencys
evDep = odep.eventDependencies
#=   println("jacobian= ",jacobian)
println("discJac= ",discJac)
println("zc_jac= ",zc_jac)
println("ZC_discJac= ",ZC_discJac)
println("evDep= ",evDep) =#
#********************************helper values*******************************
numDiscr = length(d)
states = length(initConditions)
numberZC=size(zc_jac, 1)
numEvents=numberZC*2   #each zero crossing has 1 posEv and 1 negEv
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
  push!(f, @RuntimeGeneratedFunction(odep.eqs[i].args[2])) #args[2] cuz there is extra stuff (line number in arg1)
end

zcf = Vector{Function}()
for i = 1:length(odep.zceqs)# later change to numberZC
  push!(zcf, @RuntimeGeneratedFunction(odep.zceqs[i].args[2])) #args[2] cuz there is extra stuff
end

eventf = Vector{Function}()
for i = 1:length(odep.eventEqus)# later change to numEvents
  push!(eventf, @RuntimeGeneratedFunction(odep.eventEqus[i].args[2])) 
end

#*********************************  initialize          *****************************************
t = Taylor0(order)
integratorCache=Taylor0(zeros(order),order) #for integratestate only

cacheSize=odep.cacheSize
taylorOpsCache=Array{Taylor0{Float64},1}()
for i=1:cacheSize
  push!(taylorOpsCache,Taylor0(zeros(order),order))
end

interpolationOrder = 1 # 1 means savedVars = x + t derX
#= arr = []
for i = 1:states
  push!(arr, zeros(325))
end
#savedVars = SVector{states,Array{Taylor0{Float64}}}(tuple(arr...))
savedVars = SVector{states,Array{Float64}}(tuple(arr...)) =#
savedVars = Vector{Array{Float64}}(undef, states)
for i = 1:states
  #push!(savedVars, zeros(20))# stiffness hint by                   user;& num_states
  savedVars[i]=zeros(51)  # this number can found from ft and saveat (ft/saveat=5/0.1=50) and maybe *2/3 factor (time can jump past a saving time) of user or default saveat if not provided
end

#savedTimes = Array{Float64}([initTime])
savedTimes = zeros(51) 
savedTimes[1]=initTime

  #*************************************initialize************************************
for i = 1:states# x and q get init conditions
  push!(x, Taylor0(zeros(order + 1), order) + initConditions[i])
  # push!(q, Taylor0(zeros(order), order - 1))
  push!(q, Taylor0(zeros(order+1), order))#if to use in zc function
  #push!(savedVars[i], Taylor0(zeros(interpolationOrder + 1), interpolationOrder))
  savedVars[i][1]=x[i][0]
end

for k = 1:order # compute initial derivatives for x and q (similar to a recursive way )
  for i = 1:states
    q[i].coeffs[k] = x[i].coeffs[k]  # q computed from x and it is going to be used in the next x
  end
  for i = 1:states
    clearCache(taylorOpsCache)
   
    x[i].coeffs[k+1] = ((differentiate( f[i](q,d, t + initTime,taylorOpsCache), k - 1)).coeffs[1]) / factorial(k) # /fact cuz i will store der/fac like the convention...to extract the derivatives (at endof sim) multiply by fac  derderx=coef[3]*fac(2)
   #=  f[i](integratorCache,q,d, t + initTime,taylorOpsCache)
    x[i].coeffs[k+1] = ((differentiate(integratorCache, k - 1)).coeffs[1]) / factorial(k) =#
  end
end

for i = 1:states
  tx[i] = initTime
  tq[i] = initTime
#=   for j = 1:interpolationOrder
    savedVars[i][1].coeffs[j] = x[i].coeffs[j]
  end =#
  quantum[i] = relQ * abs(x[i].coeffs[1]) 
  if quantum[i] < absQ
    quantum[i] = absQ
  end
  computeNextTime(Val(2), i, initTime, nextStateTime, x, quantum)
  clearCache(taylorOpsCache)
  computeNextInputTime(Val(2), i, initTime, 0.0,f[i](q,d,initTime + t,taylorOpsCache), nextInputTime, x,  quantum)
  #= f[i](integratorCache,q,d,initTime + t,taylorOpsCache),
  computeNextInputTime(Val(2), i, initTime, 0.0,integratorCache, nextInputTime, x,  quantum) =#

end

for i=1:numberZC
  clearCache(taylorOpsCache)
  output=zcf[i](q,d,initTime + t,taylorOpsCache) #test this evaluation
  oldsignValue[i,2]=output.coeffs[1] #value
  oldsignValue[i,1]=sign(output.coeffs[1]) #sign modify 
  computeNextEventTime(i,output,oldsignValue,initTime,  nextEventTime, quantum)
end

###################################################################################################################################################################
####################################################################################################################################################################
#---------------------------------------------------------------------------------while loop-------------------------------------------------------------------------
###################################################################################################################################################################
####################################################################################################################################################################
simt = initTime
count = 1 # not zero because intial value took 0th position
len=length(savedVars[1])

while simt < ft #&& count < 4
 
  sch = updateScheduler(nextStateTime,nextEventTime, nextInputTime)
  simt = sch[2]
  index = sch[1]
  ##########################################state########################################
  if sch[3] == :ST_STATE
    elapsed = simt - tx[index]
    integrateState(Val(2),x[index],integratorCache,elapsed)
    quantum[index] = relQ * abs(x[index].coeffs[1]) #derx=coef[2]*fac(1), derderx=coef[3]*fac(2)            
    if quantum[index] < absQ
      quantum[index] = absQ
    end
    tx[index] = simt
    for k = 1:order
      q[index].coeffs[k] = x[index].coeffs[k] #updateQ
    end
    tq[index] = simt #tq needed for higher orders
    computeNextTime(Val(2), index, simt, nextStateTime, x, quantum) #
    for i = 1:length(SD[index])
      j = SD[index][i] 
      if j != 0           
        elapsed = simt - tx[j]
        if elapsed > 0
          #"evaluate" x at new time only...derivatives get updated next using computeDerivative()
          x[j].coeffs[1] = x[j](elapsed)
          q[j].coeffs[1] = x[j].coeffs[1]
          tx[j] = simt
          tq[j] = simt
        end
        clearCache(taylorOpsCache)
         #normally and later i should update other q (integrate q=q+e derQ  for higher orders)
        integrate!(x[j],f[j](q,d,t,taylorOpsCache),x[j].coeffs[1])
        reComputeNextTime(Val(2), j, simt, nextStateTime, x, q, quantum)
      end#end if j!=0
    end#end for SD
    for i = 1:length(SZ[index])
      j = SZ[index][i] 
      if j != 0             
        #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        clearCache(taylorOpsCache)
        computeNextEventTime(j,zcf[j](q,d,simt + t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)
      end  #end if j!=0
    end#end for SZ
    ##################################input########################################
  elseif sch[3] == :ST_INPUT  # time of change has come to a state var that does not depend on anything...no one will give you a chance to change but yourself  
    elapsed = simt - tx[index]    
    clearCache(taylorOpsCache)
    computeNextInputTime(Val(2), index, simt, elapsed,  f[index](q,d,simt + t,taylorOpsCache), nextInputTime, x,  quantum)
    for i = 1:length(SD[index])
      j = SD[index][i] 
      if j != 0             
        elapsed = simt - tx[j]
        if elapsed > 0
          x[j].coeffs[1] = x[j](t + elapsed).coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
          tx[j] = simt
        end
        clearCache(taylorOpsCache)
         #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        integrate!(x[j],f[j](q,d,t,taylorOpsCache),x[j].coeffs[1]) 
        reComputeNextTime(Val(2), j, simt, nextStateTime, x, q, quantum)
      end#end if j!=0
    end#end for
    for i = 1:length(SZ[index])
      j = SZ[index][i] 
      if j != 0             
        #normally and later i should update q (integrate q=q+e derQ  for higher orders)
        clearCache(taylorOpsCache)
        computeNextEventTime(j,zcf[j](q,d,simt + t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum) 
      end  
    end
  #################################################################event########################################
  else
    #first we have a zc happened which corresponds to nexteventtime and index (one of zc) but we want also the sign in order to know ev+ or ev- 
    modifiedIndex=0
    clearCache(taylorOpsCache)
    if (zcf[index](q,d,simt + t,taylorOpsCache).coeffs[1])>0       # sign is not needed here
      modifiedIndex=2*index-1   # the  event that just occured is at  this index
    else
      modifiedIndex=2*index
    end         
    #=
      since i am checking if zc>0 then posEV Funct1 else negEvFunct2....it is better to store the whole thing as an event function
      that way i can have more complicated events without worrying about parsing...one thing: is it cheaper to store one giant event in a function
      or too many small events...i think one giant event is cheaper.
    =#
    clearCache(taylorOpsCache)
    eventf[modifiedIndex](q,d,simt + t,taylorOpsCache)
    #if a choice to use x instead of q in events, then i think there should be a q update after the eventexecuted
    nextEventTime[index]=Inf   #not tested
    for i = 1:length(HD[modifiedIndex]) # care about dependency to this event only
        j = HD[modifiedIndex][i] # 
        if j != 0
                elapsed = simt - tx[j]
                if elapsed > 0  # if event triggere by change of sign and time=now then elapsed=0
                  x[j].coeffs[1] = x[j](t + elapsed).coeffs[1] #evaluate x at new time only...derivatives get updated next using computeDerivative()
                  q[j].coeffs[1] = x[j].coeffs[1]
                  tx[j] = simt
                end
                #computeDerivative(Val(2), j,x, f[j](q(simt + t),d,simt + t)) #updating q as a vector updates all components q[0] q[1]...
                #but some q is already updated!!!causes error. if q[1] updated at 0.5 then q(t+0.5) makes another change
              clearCache(taylorOpsCache)
              integrate!(x[j],f[j](q,d,t,taylorOpsCache),x[j].coeffs[1])
              reComputeNextTime(Val(2), j, simt, nextStateTime, x, q, quantum)      
        end     
    end
    for i = 1:length(HZ[modifiedIndex])
          j = HZ[modifiedIndex][i] # this line before the next line or vice versa gave the same bench results
            if j != 0             
            #normally and later i should update q (integrate q=q+e derQ  for higher orders)
            clearCache(taylorOpsCache)
            computeNextEventTime(j,zcf[j](q,d,simt + t,taylorOpsCache),oldsignValue,simt,  nextEventTime, quantum)
          end        
    end
  end#end state input event
  if simt > savetime
      count += 1
      savetime += savetimeincrement #next savetime
      if len<count
        len=count*2
        for i=1:states
          resize!(savedVars[i],len)
        end
        resize!(savedTimes,len)
      end
      for k = 1:states
         #=  temp = Taylor0(interpolationOrder)
          for j = 1:interpolationOrder
            temp.coeffs[j] = x[k].coeffs[j]
          end
          # println("temp= ",temp)
          push!(savedVars[k], temp) =#
          #push!(savedVars[k], x[k][0])
         #=  if length(savedVars[k])<count
            resize!(savedVars[k],count*2)
          end =#
          #savedVars[k][count]=x[k][0]
         
          savedVars[k][count]=x[k][0]
    
      end
     # push!(savedTimes, simt)
     savedTimes[count]=simt
  end#end if save
end#end while
(savedTimes, savedVars)
end#end integrate