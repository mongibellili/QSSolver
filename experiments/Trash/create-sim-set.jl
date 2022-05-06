
function QSS_simGenerate(initCond::Vector{Float64},jacobian::Matrix{Float64})
    #states = computeStates(jacobian)
    states=2
    #order=getOrderfromSolverMethod(settings.solver)
    order=1
    savedVars=SVector{states,Array{Float64}}([],[])# !!!!!!!!this will have to be generated cuz of number of [][] [] []
    savedTimes=Array{Float64}([0.0])
    quantum = @MVector zeros(states)
    x = @MVector zeros(states*(order+1))#4=states*(order+1)
    q = @MVector zeros(states*(order+1))
    nextStateTime = @MVector zeros(states)
    tx = @MVector zeros(states)
    tq = @MVector zeros(states)
    minTimeValue = @MVector zeros(1)
    minIndex = MVector{1,Int}(0)
    settings = ModelSettings(initCond,jacobian,3.0,0.0,1e-4,1e-2,qss1(),saveat(0.05))
    sim= QSS_simulator(savedVars,savedTimes,quantum,x,q,tx,tq,nextStateTime,minTimeValue,minIndex,states,order)
   
    #integrate(settings,sim)
    display(sim);println()
end


function saveat(savetimeincrement::Float64)
    savetimeincrement
end 

qss1()=Val(1)
qss2()=Val(2)
qss3()=Val(3)
liqss()=Val(4)
liqss2()=Val(5)
