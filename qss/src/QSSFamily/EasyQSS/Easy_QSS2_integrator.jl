#using TimerOutputs

function QSS_integrate(::Val{2}, s::EasyQSS_data, settings::ProblemSetting,prob::easyqssProblem)
  ft = settings.finalTime
  initTime = settings.initialTime
  relQ = settings.dQrel
  absQ = settings.dQmin
  initConditions = prob.initConditions
  inputVars=prob.inputVars
  #inputJac=settings.inputJac
  solver = settings.solver
  #display(solver);println()
  jacobian = prob.jacobian
  savetimeincrement = settings.savetimeincrement
  #********************************helper values*******************************
  states = computeStates(prob.initConditions)
  #inputs = computeInputs(prob.inputVars)
  #order = getOrderfromSolverMethod(settings.solver)
  #display(order);println()
  savetime = savetimeincrement
  #dep = createDependencyMatrix(jacobian)
  #jacobian=modifyJacobian(jacobian) #either modify to transpose or create svector of svectors...to be used inside compute derivatives
  jacobian = modifyJacobian(jacobian)
  #display(jacobian);println()
  #jacobian= SMatrix{2,2,Float64}(transpose(jacobian)) #change workshop when removed. the test showed bad performance but it might be worth it for large jacobians
  dep = createDependencyMatrix(jacobian)
  #dep=((2),(1, 2))
  #dep=[[2],[1, 2]]
  #dep= @SVector[[0,2],[1,2]]
  #=   display(dep);println()
    display(typeof(dep));println() =#
  arr = []
  for i = 1:states
    push!(arr, [])
  end
  savedVars = SVector{states,Array{Float64}}(tuple(arr...))
  savedTimes = Array{Float64}([0.0])
  #*********************************data*****************************************
  quantum = s.quantum
  x = s.x
  q = s.q
  nextStateTime = s.nextStateTime
  tx = s.tx
  tq = s.tq
  #*************************************initialize************************************
  for i = 1:states
    tx[i] = initTime
    tq[i] = initTime
    x[3*i-2] = initConditions[i]
    push!(savedVars[i], initConditions[i]) #savedTimes assumed to be always initially equals zero!!!!!no
    quantum[i] = relQ * abs(x[3*i-2])
    if quantum[i] < absQ
      quantum[i] = absQ
    end
    #updateQ(solver, i, x, q, quantum)
    q[(2)*i-1] = (x[(3)*i-2])
  end

  for i = 1:states
    computeInitDerivative(i, solver, jacobian,  x, q, tx, tq,inputVars)
    #computeNextTime(solver, i, initTime, nextStateTime, x, quantum)
  end
  for i = 1:states
    #updatederQ(solver, i, x, q, quantum)
    q[(2)*i] = (x[(3)*i-1])
  end
#=   println("-----initial q------= ")
  display(q)
  println() =#
  for i = 1:states
    computeInitsecondDerivative(i, solver, jacobian,  x, q, tx, tq,inputVars)
    #computeNextTime(solver, i, initTime, nextStateTime, x, quantum)
  end
  for i = 1:states

    computeNextTime(solver, i, initTime, nextStateTime, x, quantum)
  end
#=   println("***********initial state var***********")
  display(x)
  println()
  println("***********initial next time***********")
  display(nextStateTime)
  println() =#
  t = initTime
  #**************************************integrate*************************************
  if savetimeincrement == 0
    while t < ft
      sch = updateScheduler(nextStateTime)
      t = sch[2]
      index = sch[1]
      elapsed = t - tx[index]
     # println("elapsedi= ",elapsed )
      integrateState(index, solver, elapsed, x)
      tx[index] = t
      quantum[index] = relQ * abs(x[(3)*index-2])
      if quantum[index] < absQ
        quantum[index] = absQ
      end
     # println("quantum",index,"= ",quantum[index])
      updateQ(solver, index, x, q, quantum)
     # updatederQ(solver, index, x, q, quantum)
      tq[index] = t

      computeNextTime(solver, index, t, nextStateTime, x, quantum)
      for i = 1:length(dep[index])
        if dep[index][i] != 0
          j = dep[index][i]
         # println("j= ",j)
          elapsed = t - tx[j]
          elapsedQ = t - tq[j]
         # println("elapsedj= ",elapsed )
          if elapsed > 0 #test efficiency for large number of states: cost of one extra functioncall vs cost of if statement many times
            reIntegrateState(j, solver, elapsed, elapsedQ, x, q)
            tx[j] = t
            tq[j] = t #----------------------this is needed for derx=f(q,t)   !!!!!!!!!!!!!!!!!!
          end
          computeDerivative(t,index,j, solver, jacobian,  x, q, tx, tq,inputVars)
          reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
     
        end
      end
   
    end #end while
  else
    while t < ft
      sch = updateScheduler(nextStateTime)
      t = sch[2]
      index = sch[1]
      elapsed = t - tx[index]
     # println("elapsedi= ",elapsed )
      integrateState(index, solver, elapsed, x)
      tx[index] = t
      quantum[index] = relQ * abs(x[(3)*index-2])
      if quantum[index] < absQ
        quantum[index] = absQ
      end
     # println("quantum",index,"= ",quantum[index])
      updateQ(solver, index, x, q, quantum)
     # updatederQ(solver, index, x, q, quantum)
      tq[index] = t

      computeNextTime(solver, index, t, nextStateTime, x, quantum)
      for i = 1:length(dep[index])
        if dep[index][i] != 0
          j = dep[index][i]
         # println("j= ",j)
          elapsed = t - tx[j]
          elapsedQ = t - tq[j]
         # println("elapsedj= ",elapsed )
          if elapsed > 0 #test efficiency for large number of states: cost of one extra functioncall vs cost of if statement many times
            reIntegrateState(j, solver, elapsed, elapsedQ, x, q)
            tx[j] = t
            tq[j] = t #----------------------this is needed for derx=f(q,t)   !!!!!!!!!!!!!!!!!!
          end
          computeDerivative(t,index,j, solver, jacobian,  x, q, tx, tq,inputVars)  # qss in c did this in another dependency loop
          reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
     
        end
      end

  
      if t > savetime   #if the user picks very small saveatincrements then this is inefficient; it is worse than saving original x everytime
        savetime += savetimeincrement
        for k = 1:states
          push!(savedVars[k], x[3*k-2])
        end
        push!(savedTimes, t)
      end
#=       println("***********next time***********")
      display(nextStateTime)
      println()
       println("***********state var***********")
      display(x)
      println()
      println("***********q***********")     
      display(q)
      println() 
      println("***********quantum***********")     
      display(quantum)
      println() 
      println("***********time***********")     
      display(t)
      println()  =#
    end# end while
    #print_timer()
    #=     println("***********state var***********")
        display(x);println()
        println("***********next time***********")
        display(nextStateTime);println()
          display(q);println() =#
  end# end else
  (savedTimes, savedVars)

end
