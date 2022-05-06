struct QSS_simulator{T,O,R}
    quantum :: MVector{T,Float64} 
    x :: MVector{O,Float64}
    q :: MVector{R,Float64} 
    tx ::  MVector{T,Float64} 
    tq :: MVector{T,Float64} 
    nextStateTime :: MVector{T,Float64}    
end
function QSS_simGenerate(settings :: ModelSettings)
    order=getOrderfromSolverMethod(settings.solver)
    states = computeStates(settings.initConditions)
    solver=settings.solver
    quantum = @MVector zeros(states)
    x = @MVector zeros(states*(order+1))#4=2states*(order1+1)
    q = @MVector zeros(states*(order))
    nextStateTime = @MVector zeros(states)
    tx = @MVector zeros(states)
    tq = @MVector zeros(states)
    qssSimulator= QSS_simulator(quantum,x,q,tx,tq,nextStateTime)
    QSS_integrate(solver,qssSimulator,settings)   
end
