#using TimerOutputs

function QSS_integrate(::Val{1}, s::QSS_simulator, settings::ProblemSetting,prob::qssProblem)
  #reset_timer!()
  println("normal integration! ")
  #*********************************settings*****************************************
  ft = settings.finalTime
  initTime = settings.initialTime
  relQ = settings.dQrel
  absQ = settings.dQmin
  savetimeincrement = settings.savetimeincrement
  solver = settings.solver  # if integrator going to specialize then use Val(1)
  #*********************************problem data*****************************************
  jacobian = prob.jacobian
  initConditions = prob.initConditions
  inputVars=prob.inputVars
  discreteVars=prob.discreteVars#::SVector{D,Float64} 
  zcFunctions=prob.zcFunctions#::SVector{Z,Function} 
  eventHandlerFunctions =prob.eventHandlerFunctions#::SVector{E,Function}
  SD=prob.SD#::SMatrix{T,T,Int} 
  SZ =prob.SZ#::SMatrix{T,Z,Int} # T states can influence Z zcfunctions
  HD =prob.HD#::SMatrix{D,T,Int}  # D discrete Vars can influence T derivatives
  HZ =prob.HZ#::SMatrix{D,Z,Int}  #

  #********************************helper values*******************************
  states = computeStates(initConditions)
  #inputs = computeInputs(prob.inputVars)
  events = computeEvents(eventHandlerFunctions)
  savetime = savetimeincrement
  arr = []
  for i = 1:states
    push!(arr, [])
  end
  savedVars = SVector{states,Array{Float64}}(tuple(arr...))
  savedTimes = Array{Float64}([0.0])
  numberZC=length(zcFunctions)
  oldsignValue=m1 = MMatrix{numberZC,2}(zeros(4))
  dt=1e-2
  #*********************************qss method data*****************************************
  quantum = s.quantum
  x = s.x
  q = s.q
  nextStateTime = s.nextStateTime
  tx = s.tx
  tq = s.tq
  nextEventTime = s.nextEventTime
  #*************************************initialize************************************
  for i = 1:states
    tx[i] = initTime
    tq[i] = initTime
    x[(2)*i-1] = initConditions[i]
    push!(savedVars[i], initConditions[i]) #savedTimes assumed to be always initially equals zero!!!!!no
    quantum[i] = relQ * abs(x[(2)*i-1])
    if quantum[i] < absQ
      quantum[i] = absQ
    end
    updateQ(solver, i, x, q, quantum)
  end
  #display(q);println()
  for i = 1:states
    computeDerivative(i, solver, jacobian, x, q, tx, tq,discreteVars)
    #computeNextTime(solver, i, initTime, nextStateTime, x, quantum)
  end
  for i = 1:states
    computeNextTime(solver, i, initTime, nextStateTime, x, quantum)
  end
  for i=1:events
    println("event i= ", i)
    output=zcFunctions[i](q,discreteVars,initTime,1)
    oldsignValue[i,1]=output[1] #sign
    oldsignValue[i,2]=output[2] #value
  end
  counter=MVector{1,Int}(6) # delete later ...used for limited println 
  for i=1:events
    computeNextEventTime(dt,counter,solver,i,zcFunctions,oldsignValue,dt,initTime,  nextEventTime, q,discreteVars,tq,x, quantum)
  end
 
  #**************************************integrate*************************************

  t = initTime
  if savetimeincrement == 0
    while t < ft
      sch = updateScheduler(counter,nextStateTime,nextEventTime)
      t = sch[2]
      index = sch[1]


      elapsed = t - tx[index]
      integrateState(index, solver, elapsed, x)
      tx[index] = t
      quantum[index] = relQ * abs(x[(2)*index-1])
      if quantum[index] < absQ
        quantum[index] = absQ
      end
      updateQ(solver, index, x, q, quantum)
      tq[index] = t
      computeNextTime(solver, index, t, nextStateTime, x, quantum)
      for i = 1:size(SD,2)  # n columns
        if SD[index,i] != 0
          j = SD[index,i]
          elapsed = t - tx[j]
          if elapsed > 0 #test efficiency for large number of states: cost of one extra functioncall vs cost of if statement many times
            integrateState(j, solver, elapsed, x)
            tx[j] = t
                  for i=1:size(SZ,2)  # this code is added here cuz x2 became so fast that x1 only change here (nexttime of 2 is always the min)
                    #println("index qq SZ[]= ",index)  index stuck at 2
                    if SZ[j,i] != 0
                    # println("index SZ= ",index)
                      j = SZ[j,i] # zc function number j depends on x[index]
                      #normally and later i should update q (integrate q=q+e derQ  for higher orders)
                      computeNextEventTime(elapsed,counter,solver,j,zcFunctions,oldsignValue,dt,t,  nextEventTime, q,discreteVars,tq,x, quantum)
                      #println("there is a state var that influenced zc ")
                    end  
                  end
          end
          computeDerivative(j, solver, jacobian, x, q, tx, tq,discreteVars)
          reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
        end
      end

    end #end while
  else
    while t < ft
      sch = updateScheduler(counter,nextStateTime,nextEventTime)
      t = sch[2]
      index = sch[1]
      if sch[3] ==1
        #decide here on type state or event
        # if state ----------------------------------------------------------

        elapsed = t - tx[index]
        integrateState(index, solver, elapsed, x)
        tx[index] = t
        quantum[index] = relQ * abs(x[(2)*index-1])
        if quantum[index] < absQ
          quantum[index] = absQ
        end
        updateQ(solver, index, x, q, quantum)
        tq[index] = t
        computeNextTime(solver, index, t, nextStateTime, x, quantum)
        for i = 1:size(SD,2)  # n columns
          if SD[index,i] != 0
            j = SD[index,i]
            elapsed = t - tx[j]
            if elapsed > 0
              integrateState(j, solver, elapsed, x)
              tx[j] = t
              #= for i=1:size(SZ,2)  # this code is added here cuz x2 became so fast that x1 only change here (nexttime of 2 is always the min)
                #println("index qq SZ[]= ",index)  index stuck at 2
                if SZ[j,i] != 0
                # println("index SZ= ",index)
                  j = SZ[j,i] # zc function number j depends on x[index]
                  #normally and later i should update q (integrate q=q+e derQ  for higher orders)
                  computeNextEventTime(elapsed,counter,solver,j,zcFunctions,oldsignValue,dt,t,  nextEventTime, q,discreteVars,tq,x, quantum)
                  #println("there is a state var that influenced zc ")
                end  
              end =#




            end
            computeDerivative(j, solver, jacobian, x, q, tx, tq,discreteVars)
            reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
          end
        end
        for i=1:size(SZ,2)
          #println("index qq SZ[]= ",index)  index stuck at 2
          if SZ[index,i] != 0
           # println("index SZ= ",index)
            j = SZ[index,i] # zc function number j depends on x[index]
            #normally and later i should update q (integrate q=q+e derQ  for higher orders)
            computeNextEventTime(elapsed,counter,solver,j,zcFunctions,oldsignValue,dt,t,  nextEventTime, q,discreteVars,tq,x, quantum)
            #println("there is a state var that influenced zc ")
          end  
        end
      
        #---------------------end of state
        #if event
    else
      #if sch[3]==2
        println("handling event ",index," at time= ",t)
        eventHandlerFunctions[index](q,discreteVars,tq)
        nextEventTime[index]=Inf
        for i = 1:size(HD,1)  # n rows
              for k=1:size(HD,2)
                if HD[i,k] != 0
                    j = HD[i,k]
                    elapsed = t - tx[j]
                    if elapsed > 0
                      integrateState(j, solver, elapsed, x)
                      tx[j] = t
                    end
                    computeDerivative(j, solver, jacobian, x, q, tx, tq,discreteVars)
                    reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
                end
              end
        end
        for i=1:size(HZ,1)
            for k=1:size(HZ,2)
              if HZ[i,k] != 0
                j = HZ[i,k] # zc function number j depends on x[index]
                #normally and later i should update q (integrate q=q+e derQ  for higher orders)
                computeNextEventTime(dt,counter,solver,j,zcFunctions,oldsignValue,dt,currentTime,  nextEventTime, q,discreteVars,tq,x, quantum)
                println("there is a discrete var that influenced zc ")
              end  
            end
        end





      end
        #------------------end of event

      if t > savetime   #if the user picks very small saveatincrements then this is inefficient; it is worse than saving original x everytime
        savetime += savetimeincrement
        for k = 1:states
          push!(savedVars[k], x[(2)*k-1])
        end
        push!(savedTimes, t)
      end
            #= println("***********state var***********")
      display(x)
      println() =#
      #=    println("***********next time***********")
      display(nextStateTime) =#
      #= println()
      display(q)
      println() =#
    end# end while
  end# end else
  (savedTimes, savedVars)
  #print_timer()
end
