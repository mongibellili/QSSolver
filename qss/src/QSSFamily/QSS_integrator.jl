#using TimerOutputs

function QSS_integrate(s::QSS_simulator,settings::ModelSettings)
  #reset_timer!()
  #*********************************settings*****************************************
  ft = settings.finalTime
  initTime = settings.initialTime
  relQ = settings.dQrel
  absQ = settings.dQmin
  initConditions = settings.initConditions
  solver = settings.solver
  jacobian = settings.jacobian
  savetimeincrement=settings.savetimeincrement
  #********************************helper values*******************************
  states = computeStates(settings.initConditions) 
  order=getOrderfromSolverMethod(settings.solver)
  savetime=savetimeincrement
  #dep = createDependencyMatrix(jacobian)
  jacobian=modifyJacobian(jacobian) #either modify to transpose or create svector of svectors...to be used inside compute derivatives
  #jacobian= SMatrix{2,2,Float64}(transpose(jacobian)) #change workshop when removed. the test showed bad performance but it might be worth it for large jacobians
  dep = createDependencyMatrix(jacobian)
  #dep=((2),(1, 2))
  #dep=[[2],[1, 2]]
  #dep= @SVector[[2],[1,2]]
#=   display(dep);println()
  display(typeof(dep));println() =#
  arr=[]
  for i = 1:states 
      push!(arr,[])        
  end
  savedVars=SVector{states,Array{Float64}}(tuple(arr...))
  savedTimes=Array{Float64}([0.0])
  #*********************************data*****************************************
  quantum=s.quantum 
  x=s.x 
  q=s.q 
  nextStateTime=s.nextStateTime 
  tx=s.tx
  tq=s.tq 
  #*************************************initialize************************************
  for i = 1:states
    tx[i] = initTime
    tq[i] = initTime
    x[(order+1)*i-order] = initConditions[i]
    push!(savedVars[i],initConditions[i]) #savedTimes assumed to be always initially equals zero!!!!!no
    quantum[i] = relQ * abs(x[(order+1)*i-order])
    if quantum[i] < absQ
      quantum[i] = absQ
    end
    updateQ(solver, i, x, q, quantum)
  end
  for i = 1:states
    computeDerivative( i, order, jacobian, x, q, tx, tq)
    computeNextTime(solver, i, initTime, nextStateTime, x, quantum)
  end
  sch=updateScheduler(nextStateTime)
  t = sch[2]
  index = sch[1]
  #**************************************integrate*************************************
  if savetimeincrement ==0
    while t < ft
      elapsed = t - tx[index]
      integrateState(index, solver, elapsed, x)
       tx[index] = t
      quantum[index] = relQ * abs(x[(order+1)*index-order])
      if quantum[index] < absQ
        quantum[index] = absQ
      end
       updateQ(solver, index, x, q, quantum)
      tq[index] = t
       computeNextTime(solver, index, t, nextStateTime, x, quantum)
      for i = 1:length(dep[index])
         if dep[index][i]!=0
            j = dep[index][i]
         
        elapsed = t - tx[j]
          integrateState(j, solver, elapsed, x)
         tx[j] = t
         computeDerivative( j, order, jacobian, x, q, tx, tq)
         reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
         end
      end
        sch=updateScheduler(nextStateTime)
      t = sch[2]
      index = sch[1]
    end
  else
    while t < ft
      elapsed = t - tx[index]
      integrateState(index, solver, elapsed, x)
       tx[index] = t
      quantum[index] = relQ * abs(x[(order+1)*index-order])
      if quantum[index] < absQ
        quantum[index] = absQ
      end
       updateQ(solver, index, x, q, quantum)
      tq[index] = t
       computeNextTime(solver, index, t, nextStateTime, x, quantum)
      for i = 1:length(dep[index])
         if dep[index][i]!=0
            j = dep[index][i]
         
        elapsed = t - tx[j]
          integrateState(j, solver, elapsed, x)
         tx[j] = t
         computeDerivative( j, order, jacobian, x, q, tx, tq)
         reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
         end
      end
        sch=updateScheduler(nextStateTime)
      t = sch[2]
      index = sch[1]
      if t>savetime   #if the user picks very small saveatincrements then this is inefficient; it is worse than saving original x everytime
        savetime+=savetimeincrement
        for k=1:states
          push!(savedVars[k],x[(order+1)*k-order])
        end
        push!(savedTimes,t)
      end
    end
  end
 (savedTimes,savedVars)
 #print_timer()
end
