struct QSS_simulator{T,O}
#=     savedVars::SVector{T,Array{Float64}}
    savedTimes::Array{Float64} =#
    quantum :: MVector{T,Float64} 
    x :: MVector{O,Float64}
    q :: MVector{O,Float64} 
    tx ::  MVector{T,Float64} 
    tq :: MVector{T,Float64} 
    nextStateTime :: MVector{T,Float64}    
end
function QSS_simGenerate(settings :: ModelSettings)
    order=getOrderfromSolverMethod(settings.solver)
    states = computeStates(settings.initConditions)
#=     arr=[]
    for i = 1:states 
        push!(arr,[])        
    end
    savedVars=SVector{states,Array{Float64}}(tuple(arr...))
    savedTimes=Array{Float64}([0.0]) =#
    quantum = @MVector zeros(states)
    x = @MVector zeros(states*(order+1))#4=states*(order+1)
    q = @MVector zeros(states*(order+1))
    nextStateTime = @MVector zeros(states)
    tx = @MVector zeros(states)
    tq = @MVector zeros(states)
    qssSimulator= QSS_simulator(quantum,x,q,tx,tq,nextStateTime)
    QSS_integrate(qssSimulator,settings)   
end
