
  function QSS_integrate(qssSimulator:: QSS_simulator)
   
    #scheduler=file that contains a function to find min and update time
    #quantizer=file that contains 3 functions: computeNext, recomputeNext, updateQ
    #framework= file that contains 3 functions: integrateState, compute 1 derivative, compute dependent derivatives
   
    #----------setting------
    order=qssSimulator.settings.order
    coeffx=order+1
    ft=qssSimulator.settings.finalTime
    initTime=qssSimulator.settings.initialTime
    relQ=qssSimulator.settings.dQrel
    absQ=qssSimulator.settings.dQmin
    initConditions=qssSimulator.settings.initConditions

    #---------data---------
    qssdata=qssSimulator.qssData
    states=qssdata.states 
    x=qssdata.x
    q=qssdata.q
    quantum=qssdata.quantum
    #--------time--------
    qsstime=qssSimulator.qssTime
    t=qsstime.time # make sure to update time after change...passbyValue
   # println("initial time= ", t)
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
        #push!(q[(order+1)*i-order],initConditions[i]) #normally q should be updated through quantizer
          q[(order+1)*i-order]=initConditions[i]
    #--compute initial deltaQ  
      quantum[i]=relQ * abs(x[(order+1)*i-order]) 
        if quantum[i] < absQ
              quantum[i]=absQ
        end
        
    end
    
    #-----------initial derivatives: ask framework to compute derivatives=f(x,t)
    #----------initial nextTimes: ask quantizer to compute nextChangeTime
    for i = 1:states
        computeDerivative(states,i,order,qssmodel,x,q,tx,tq)
        computeNextTime(quantizer,i,t,tn,x,quantum)
    end

   
    #----------ask scheduler to update time and finds next minTime and minIndex
    updateScheduler(qsstime)
    t=qsstime.time
  
    index=qsstime.minIndex
    #**************************************integrate*************************************

    while t < ft
       
       # println("index $index at time $t ")
        elapsed=t-tx[index]
        integrateState(index,order,elapsed,x)
          tx[index]=t
           quantum[index]=relQ * abs( x[(order+1)*index-order] ) 
        if quantum[index] < absQ
              quantum[index]=absQ
        end
        updateQ(quantizer,index,x,q,quantum)# the whole quantum is not needed
       # if length(tq[index])!=0
         #   pop!(tq[index])
       # end
        #push!(tq[index],t)
         tq[index]=t
        computeNextTime(quantizer,index,t,tn,x,quantum)
       
        
        for i=1:length(dep[index])
             j=qssmodel.dep[index][i]
                elapsed=t-tx[j]
            integrateState(j,order,elapsed,x)
           # println("derx= ", x[(order+1)*j])
            #if j != index
            
                tx[j]=t
            #end
            computeDerivative(states,j,order,qssmodel,x,q,tx,tq)
            reComputeNextTime(quantizer,j,t,tn,x,q,quantum)
        end
        #computeDerivative(states,index,order,qssmodel,x,q,tx,tq)
        #reComputeNextTime(quantizer,index,t,tn,x,q,quantum)
        #println("x= ", x[(order+1)*index-order])
       # println("derx= ", x[(order+1)*index])
        updateScheduler(qsstime)
        t=qsstime.time
        index=qsstime.minIndex
        


    end
  
   
    
            
   #= 
    println("derx1= ", x[(order+1)*1])
    println("derx2= ", x[(order+1)*2])
    println("q1= ", q[1])
    println("q2= ", q[3])
    println("tq1= ", tq[1])
    println("tq2= ", tq[2])
    =#
   # println("tx1= ", tx[1])
   # println("x1= ", x[1])
   

end


# where did i stop: create func integrate state