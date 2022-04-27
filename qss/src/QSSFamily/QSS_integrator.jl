
function createDependencyMatrix(jac :: SMatrix{2,2,Float64}  )
  # extract dependency matrix from jac
 #=  epselon=1e-6
  nRows=size(jac,1)
  nColumns=size(jac,2)
  dep=MVector{nColumns,Array{Int}}([],[])###################################optimize meeeeeeee
  for j=1:nColumns
       #dep[j]=Array{Float64}[]
       for i=1:nRows
           if jac[i,j] < -epselon ||  jac[i,j] > epselon # different than zero
               push!(dep[j],i)
           end
       end
   end
  # dep=@SVector{nColumns,}
  return dep   =#
  epselon=1e-6
  nRows=size(jac,1)
  nColumns=size(jac,2)
  dep=Vector{Array{Int}}(undef, nColumns)

  for j=1:nColumns
       dep[j]=Array{Float64}[]
       for i=1:nRows
           if jac[i,j] < -epselon ||  jac[i,j] > epselon # different than zero
               push!(dep[j],i)
           end
       end
   end
   

   return dep
end
function computeStates(jacobian :: SMatrix{2,2,Float64})
   # return number of rows
   #return size(jacobian,1)
   return 2
end
function getOrderfromSolverMethod(::Val{1})
    1
end
function modifyJacobian(jacobian :: SMatrix{2,2,Float64})
  #test use of svector of svectos
end
function QSS_integrate(settings::ModelSettings)
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
  order=getOrderfromSolverMethod(solver)
  savetime=savetimeincrement
  dep = createDependencyMatrix(jacobian)
  states = computeStates(jacobian)
  #modifiedJac=modifyJacobian(jacobian) #either modify to transpose or create svector of svectors...to be used inside compute derivatives
  jacobian= SMatrix{2,2,Float64}(transpose(jacobian))
  savedVars=SVector{2,Array{Float64}}([],[])# !!!!!!!!this will have to be generated cuz of number of [][] [] []
  savedTimes=Array{Float64}([0.0])
  #*********************************data*****************************************
  quantum = @MVector zeros(2)
  x = @MVector zeros(4)
  q = @MVector zeros(4)
  nextStateTime = @MVector zeros(2)
  tx = @MVector zeros(2)
  tq = @MVector zeros(2)
  #reset_timer!() 
  minTimeValue = @MVector zeros(1)
  minTimeValue[1] = initTime
  minIndex = MVector{1,Int}(0) #minindex = 0 not tested meaning of =0
  #*************************************initialize************************************
  for i = 1:states
    tx[i] = initTime
    tq[i] = initTime
    x[(order+1)*i-order] = initConditions[i]
    push!(savedVars[i],initConditions[i])
    quantum[i] = relQ * abs(x[(order+1)*i-order])
    if quantum[i] < absQ
      quantum[i] = absQ
    end
    updateQ(solver, i, x, q, quantum)
  end
  for i = 1:states
    computeDerivative(states, i, order, jacobian, x, q, tx, tq)
    computeNextTime(solver, i, minTimeValue[1], nextStateTime, x, quantum)
  end
  updateScheduler(states, nextStateTime, minTimeValue, minIndex)
  t = minTimeValue[1]
  index = minIndex[1]
  #**************************************integrate*************************************
  if savetimeincrement ==0
    while t < ft
      elapsed = t - tx[index]
      integrateState(index, order, elapsed, x)
      tx[index] = t
      quantum[index] = relQ * abs(x[(order+1)*index-order])
      if quantum[index] < absQ
        quantum[index] = absQ
      end
      updateQ(solver, index, x, q, quantum)
      tq[index] = t
      computeNextTime(solver, index, t, nextStateTime, x, quantum)
      for i = 1:length(dep[index])
        j = dep[index][i]
        elapsed = t - tx[j]
        integrateState(j, order, elapsed, x)
        tx[j] = t
        computeDerivative(states, j, order, jacobian, x, q, tx, tq)
        reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
      end
      updateScheduler(states, nextStateTime, minTimeValue, minIndex)
      t = minTimeValue[1]
      index = minIndex[1]
    end
  else
    while t < ft
      elapsed = t - tx[index]
      integrateState(index, order, elapsed, x)
      tx[index] = t
      quantum[index] = relQ * abs(x[(order+1)*index-order])
      if quantum[index] < absQ
        quantum[index] = absQ
      end
      updateQ(solver, index, x, q, quantum)
      tq[index] = t
      computeNextTime(solver, index, t, nextStateTime, x, quantum)
      for i = 1:length(dep[index])
        j = dep[index][i]
        elapsed = t - tx[j]
        integrateState(j, order, elapsed, x)
        tx[j] = t
        computeDerivative(states, j, order, jacobian, x, q, tx, tq)
        reComputeNextTime(solver, j, t, nextStateTime, x, q, quantum)
      end
      updateScheduler(states, nextStateTime, minTimeValue, minIndex)
      t = minTimeValue[1]
      index = minIndex[1]
      if t>savetime
        savetime+=savetimeincrement
        for k=1:states
          push!(savedVars[k],x[(order+1)*k-order])
          
        end
        push!(savedTimes,t)
      end
    end
  end
 (savedTimes,savedVars)
end
