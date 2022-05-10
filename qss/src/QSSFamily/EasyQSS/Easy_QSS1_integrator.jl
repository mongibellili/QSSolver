#using TimerOutputs

function QSS_integrate(::Val{1}, s::EasyQSS_simulator, settings::ProblemSetting,prob::easyqssProblem)
  println("easy integration!")
  #reset_timer!()
  #*********************************settings*****************************************
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
  #nextInputTime = s.nextInputTime
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
    computeDerivative(i, solver, jacobian, x, q, tx, tq,inputVars)
    #computeNextTime(solver, i, initTime, nextStateTime, x, quantum)
  end
  for i = 1:states

    computeNextTime(solver, i, initTime, nextStateTime, x, quantum)
  end
#=   for i = 1:inputs

    computeNextInput(solver, i, initTime, nextStateTime, x, quantum)
  end =#
  #=   println("***********state var***********")
    display(x);println()
    println("***********next time***********")
    display(nextStateTime);println() =#
  #=   sch=updateScheduler(nextStateTime)
    t = sch[2]
    index = sch[1] =#
  #**************************************integrate*************************************
  t = initTime
  if savetimeincrement == 0
    while t < ft
      sch = updateScheduler(nextStateTime)
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
      for i = 1:length(dep[index])
        if dep[index][i] != 0
          j = dep[index][i]

          elapsed = t - tx[j]
          if elapsed > 0 #test efficiency for large number of states: cost of one extra functioncall vs cost of if statement many times
            integrateState(j, solver, elapsed, x)
            tx[j] = t
          end
          computeDerivative(j, solver, jacobian, x, q, tx, tq,inputVars)
          reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
        end
      end

    end
  else
    while t < ft
      sch = updateScheduler(nextStateTime)
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
      for i = 1:length(dep[index])
        if dep[index][i] != 0
          j = dep[index][i]
          elapsed = t - tx[j]
          if elapsed > 0
            integrateState(j, solver, elapsed, x)
            tx[j] = t
          end
          computeDerivative(j, solver, jacobian, x, q, tx, tq,inputVars)
          reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
        end
      end

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
