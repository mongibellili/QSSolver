#= using TimerOutputs
using InteractiveUtils =#
function QSS_integrate(::Val{1}, s::MediumQSS_data, settings::ProblemSetting, prob::NoJacMediumProblem)
  #reset_timer!()
  #*********************************settings*****************************************
  ft = settings.finalTime
  initTime = settings.initialTime
  relQ = settings.dQrel
  absQ = settings.dQmin
  initConditions = prob.initConditions
  solver = settings.solver
  #jacobian = prob.jacobian
  #display(jacobian[1])
  #make the functions available?
  #for i=1:length(jacobian)
    
 
  #end
  SD = prob.SD
  dep = modifyJacobian2(SD)
  savetimeincrement = settings.savetimeincrement
  #********************************helper values*******************************
  states = computeStates(prob.initConditions)
  savetime = savetimeincrement
  
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
    x[(2)*i-1] = initConditions[i]
    push!(savedVars[i], initConditions[i]) #savedTimes assumed to be always initially equals zero!!!!!no
    quantum[i] = relQ * abs(x[(2)*i-1])
    if quantum[i] < absQ
      quantum[i] = absQ
    end
    updateQ(solver, i, x, q, quantum)
  end
 #=  f(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int) where {R,T}=q[2] 
  function computeDerivative(index::Int, ::Val{1}   ,x::MVector{O,Float64} , q::MVector{T,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64},   jacobian::SVector{T,Function} ) where{T,O} 
    #f(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int) where {R,T}=jacobian[index](q,tq,1)
     
    computeDerivative(index, Val(1)  ,x , q, tx, tq,f ) 
  end =#
  for i = 1:states
    computeDerivative(i, solver, x, q, tx, tq)
  end
  for i = 1:states
    computeNextTime(solver, i, initTime, nextStateTime, x, quantum)
  end
  
  #**************************************integrate*************************************
  t = initTime
  tem=0
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
      #=       for i = 1:size(SD,2)  # n columns
              if SD[index,i] != 0
                j = SD[index,i] =#
      for i = 1:length(dep[index])
        if dep[index][i] != 0
          j = dep[index][i]
          elapsed = t - tx[j]
          if elapsed > 0 #test efficiency for large number of states: cost of one extra functioncall vs cost of if statement many times
            integrateState(j, solver, elapsed, x)
            tx[j] = t
          end
         computeDerivative(j, solver, x, q, tx, tq)
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
      #@timeit "compute next" 
      computeNextTime(solver, index, t, nextStateTime, x, quantum)
      #=       for i = 1:size(SD,2)  # n columns
              if SD[index,i] != 0
                j = SD[index,i] =#
      for i = 1:length(dep[index])
        j = dep[index][i]
        if j != 0

          elapsed = t - tx[j]
          if elapsed > 0
            #@timeit "integrate state" 
            integrateState(j, solver, elapsed, x)
            tx[j] = t
          end
          #gh=tem
          #@timeit "compu deri" 
          computeDerivative(j, solver, x, q, tx, tq)
          #x[(2)*j]=tem
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
    end# end while
  end# end else
  #print_timer()
  (savedTimes, savedVars)
  
end


 function funcX(::Val{1} ,q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int)where{R,T}
   q[2]
end
function funcX(::Val{2} ,q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int)where{R,T} 
  -q[1]-q[2]
end

function funcX1(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int)where{R,T}
  q[2]
end
function funcX2(q::MVector{R,Float64} ,tq::MVector{T,Float64} ,order::Int)where{R,T} 
 -q[1]-q[2]
end