#hold helper datastructures needed for simulation
struct QSS_data{T,Z}
    quantum :: Vector{Float64} 
 #=    x :: MVector{O,Float64}
    q :: MVector{R,Float64}  =#
    x :: Vector{Taylor0{Float64}}
    derx :: Vector{Taylor0{Float64}}
    q :: Vector{Taylor0{Float64}}
    tx ::  MVector{T,Float64} 
    tq :: MVector{T,Float64} 
    nextStateTime :: MVector{T,Float64}    
    nextInputTime :: Vector{Float64}  
    nextEventTime :: MVector{Z,Float64}  
    order::Int
end


function QSS_Solve(settings :: SimSetting,prob::NLODEProblem)
    solver=settings.solver
    order=getOrderfromSolverMethod(solver)
    states = computeStates(prob.initConditions)
   # ZCevents = computeStates(prob.ZCinputVars)
   #ZCevents=0
   numberZC=size(prob.ZC_jacobian, 1)  #later test prob.Z
    
    quantum =  zeros(states)
#=     x = @MVector zeros(states*(order+1))#4=2states*(order1+1)
    q = @MVector zeros(states*(order)) =#
    x = Vector{Taylor0{Float64}}()
    derx = Vector{Taylor0{Float64}}()
    q = Vector{Taylor0{Float64}}()
    nextStateTime = @MVector zeros(states)
    nextInputTime =  zeros(states)
    tx = @MVector zeros(states)
    tq = @MVector zeros(states)
    nextEventTime=@MVector zeros(numberZC)
    qssdata= QSS_data(quantum,x,derx,q,tx,tq,nextStateTime,nextInputTime ,nextEventTime,order)
    QSS_integrate(solver,qssdata,settings,prob)   # solver=Val(\d) needed to specialize the integrator
    #return nothing be careful to add this
end

