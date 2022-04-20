
#=
#integrate without saving
  #using TimerOutputs
  function QSS_integrate(qssSimulator:: QSS_simulator)
  #  reset_timer!()
    #scheduler=file that contains a function to find min and update time
    #quantizer=file that contains 3 functions: computeNext, recomputeNext, updateQ
    #framework= file that contains 3 functions: integrateState, compute 1 derivative, compute dependent derivatives

    #----------setting------
    order=qssSimulator.settings.order
    casesVect = Vector{Float64}(undef, 2) #to be used for cases when recomputing nexttime
    ft=qssSimulator.settings.finalTime
    initTime=qssSimulator.settings.initialTime
    relQ=qssSimulator.settings.dQrel
    absQ=qssSimulator.settings.dQmin
    initConditions=qssSimulator.settings.initConditions
    solver=qssSimulator.settings.solver

    #---------data---------
    qssdata=qssSimulator.qssData
    states=qssdata.states 
    x=qssdata.x
    q=qssdata.q
    quantum=qssdata.quantum
    #--------time--------
    qsstime=qssSimulator.qssTime
    t=qsstime.time # make sure to update time after change...passbyValue

    tx=qsstime.tx
    tq=qsstime.tq
    tn=qsstime.nextStateTime
    #---------model-------
    qssmodel=qssSimulator.qssModel
    dep=qssmodel.dep
   # display(qssmodel.jacobian)
   quantizer=Quantizer(qssSimulator)
    #*************************************initialize************************************
    for i = 1:states
    #---initial time for each variable
         # tx[i]=Array{Float64}[]
        tx[i]=initTime
        #tq[i]=Array{Float64}[]
        #push!(tq[i],initTime)
           tq[i]=initTime

    # initial condition for each variable
        #(order+1)*i-order means only the state, no derivatives
        x[(order+1)*i-order]=initConditions[i]#push===fill in first element of arr     
          #--compute initial deltaQ  
      quantum[i]=relQ * abs(x[(order+1)*i-order]) 
      if quantum[i] < absQ
            quantum[i]=absQ
      end      
        #push!(q[(order+1)*i-order],initConditions[i]) #normally q should be updated through quantizer
          #q[(order+1)*i-order]=initConditions[i]
          updateQ(solver,quantizer,i,x,q,quantum)


    end

    #-----------initial derivatives: ask framework to compute derivatives=f(x,t)
    #----------initial nextTimes: ask quantizer to compute nextChangeTime
    for i = 1:states
        computeDerivative(states,i,order,qssmodel,x,q,tx,tq)
         computeNextTime(solver,quantizer,i,t,tn,x,quantum)
        #reComputeNextTime(quantizer,i,t,tn,x,q,quantum,casesVect)
    end


    #----------ask scheduler to update time and finds next minTime and minIndex
    updateScheduler(qsstime)
    t=qsstime.time

    index=qsstime.minIndex
    #**************************************integrate*************************************

   while t < ft

       # println("index $index at time $t ")
        elapsed=t-tx[index]
        #@timeit "integrate state" integrateState(index,order,elapsed,x)
         integrateState(index,order,elapsed,x)
          tx[index]=t
           quantum[index]=relQ * abs( x[(order+1)*index-order] ) 
        if quantum[index] < absQ
              quantum[index]=absQ
        end
       # @timeit "updateQ" updateQ(quantizer,index,x,q,quantum)# the whole quantum is not needed
        updateQ(solver,quantizer,index,x,q,quantum)
       # if length(tq[index])!=0
         #   pop!(tq[index])
       # end
        #push!(tq[index],t)
         tq[index]=t
        # @timeit "compute next" computeNextTime(quantizer,index,t,tn,x,quantum)
         computeNextTime(solver,quantizer,index,t,tn,x,quantum)


        for i=1:length(dep[index])
             j=qssmodel.dep[index][i]
                elapsed=t-tx[j]
           integrateState(j,order,elapsed,x)
           # println("derx= ", x[(order+1)*j])
            #if j != index

                tx[j]=t
            #end
          #  @timeit "compute deriv"  computeDerivative(states,j,order,qssmodel,x,q,tx,tq)
           computeDerivative(states,j,order,qssmodel,x,q,tx,tq)
           #@timeit "recompute next" reComputeNextTime(quantizer,j,t,tn,x,q,quantum)

            reComputeNextTime(solver,quantizer,j,t,tn,x,q,quantum,casesVect)
        end
        #computeDerivative(states,index,order,qssmodel,x,q,tx,tq)
        #reComputeNextTime(quantizer,index,t,tn,x,q,quantum)
        #println("x= ", x[(order+1)*index-order])
       # println("derx= ", x[(order+1)*index])
        updateScheduler(qsstime)

        t=qsstime.time
        index=qsstime.minIndex
       # print_timer()


   end




  #=
    println("derx1= ", x[(order+1)*1])
    println("derx2= ", x[(order+1)*2])
    println("q1= ", q[1])
    println("q2= ", q[3])
    println("tq1= ", tq[1])
    println("tq2= ", tq[2])

    println("tx1= ", tx[1])
    println("x1= ", x[1])
    println("tx2= ", tx[2])
    println("x2= ", x[2])
    println("tn1= ", tn[1])
    println("tn2= ", tn[2])
  =#
 # println("qssdata= ", qssdata)
 # println("qsstime= ", qsstime)
 # println("qssmodel= ", qssmodel)

end



=#

#using TimerOutputs
function QSS_integrate(qssSimulator::QSS_simulator)
  #  reset_timer!()
  #scheduler=file that contains a function to find min and update time
  #quantizer=file that contains 3 functions: computeNext, recomputeNext, updateQ
  #framework= file that contains 3 functions: integrateState, compute 1 derivative, compute dependent derivatives

  #----------setting------
  order = qssSimulator.settings.order
  casesVect = Vector{Float64}(undef, 2) #to be used for cases when recomputing nexttime
  ft = qssSimulator.settings.finalTime
  initTime = qssSimulator.settings.initialTime
  relQ = qssSimulator.settings.dQrel
  absQ = qssSimulator.settings.dQmin
  initConditions = qssSimulator.settings.initConditions
  solver = qssSimulator.settings.solver

  #---------data---------
  qssdata = qssSimulator.qssData
  states = qssdata.states
  x = qssdata.x
  q = qssdata.q
  quantum = qssdata.quantum
  #--------time--------
  qsstime = qssSimulator.qssTime
  t = qsstime.time # make sure to update time after change...passbyValue

  tx = qsstime.tx
  tq = qsstime.tq
  tn = qsstime.nextStateTime
  #---------model-------
  qssmodel = qssSimulator.qssModel
  dep = qssmodel.dep
  # display(qssmodel.jacobian)
  quantizer = Quantizer(qssSimulator)
  #*************************************initialize************************************
  for i = 1:states
    #---initial time for each variable
    tx[i] = Array{Float64}[]

    #tq[i]=Array{Float64}[]
    push!(tx[i], initTime)
    tq[i] = initTime

    # initial condition for each variable
    #(order+1)*i-order means only the state, no derivatives
    push!(x[(order+1)*i-order], initConditions[i])#push===fill in first element of arr     
    #--compute initial deltaQ  
    quantum[i] = relQ * abs(last(x[(order+1)*i-order]))
    if quantum[i] < absQ
      quantum[i] = absQ
    end
    #push!(q[(order+1)*i-order],initConditions[i]) #normally q should be updated through quantizer
    #q[(order+1)*i-order]=initConditions[i]
    updateQ(solver, quantizer, i, x, q, quantum)


  end

  #-----------initial derivatives: ask framework to compute derivatives=f(x,t)
  #----------initial nextTimes: ask quantizer to compute nextChangeTime
  for i = 1:states
    computeDerivative(states, i, order, qssmodel, x, q, tx, tq)
    computeNextTime(solver, quantizer, i, t, tn, x, quantum)
    #reComputeNextTime(quantizer,i,t,tn,x,q,quantum,casesVect)
  end


  #----------ask scheduler to update time and finds next minTime and minIndex
  updateScheduler(qsstime)
  t = qsstime.time

  index = qsstime.minIndex
  #**************************************integrate*************************************

  while t < ft

    # println("index $index at time $t ")
    elapsed = t - last(tx[index])
    #@timeit "integrate state" integrateState(index,order,elapsed,x)
    integrateState(index, order, elapsed, x)
    push!(tx[index], t)
    quantum[index] = relQ * abs(last(x[(order+1)*index-order]))
    if quantum[index] < absQ
      quantum[index] = absQ
    end
    # @timeit "updateQ" updateQ(quantizer,index,x,q,quantum)# the whole quantum is not needed
    updateQ(solver, quantizer, index, x, q, quantum)
    # if length(tq[index])!=0
    #   pop!(tq[index])
    # end
    #push!(tq[index],t)
    tq[index] = t
    # @timeit "compute next" computeNextTime(quantizer,index,t,tn,x,quantum)
    computeNextTime(solver, quantizer, index, t, tn, x, quantum)


    for i = 1:length(dep[index])
      j = qssmodel.dep[index][i]
      elapsed = t - last(tx[j])
      integrateState(j, order, elapsed, x)
      # println("derx= ", x[(order+1)*j])
      #if j != index

      push!(tx[j], t)
      #end
      #  @timeit "compute deriv"  computeDerivative(states,j,order,qssmodel,x,q,tx,tq)
      computeDerivative(states, j, order, qssmodel, x, q, tx, tq)
      #@timeit "recompute next" reComputeNextTime(quantizer,j,t,tn,x,q,quantum)

      reComputeNextTime(solver, quantizer, j, t, tn, x, q, quantum, casesVect)
    end
    #computeDerivative(states,index,order,qssmodel,x,q,tx,tq)
    #reComputeNextTime(quantizer,index,t,tn,x,q,quantum)
    #println("x= ", x[(order+1)*index-order])
    # println("derx= ", x[(order+1)*index])
    updateScheduler(qsstime)

    t = qsstime.time
    index = qsstime.minIndex
    # print_timer()


  end




  #=
    println("derx1= ", x[(order+1)*1])
    println("derx2= ", x[(order+1)*2])
    println("q1= ", q[1])
    println("q2= ", q[3])
    println("tq1= ", tq[1])
    println("tq2= ", tq[2])

    println("tx1= ", tx[1])
    println("x1= ", x[1])
    println("tx2= ", tx[2])
    println("x2= ", x[2])
    println("tn1= ", tn[1])
    println("tn2= ", tn[2])
  =#
  # println("qssdata= ", qssdata)
  # println("qsstime= ", qsstime)
  # println("qssmodel= ", qssmodel)

end


# where did i stop: create func integrate state