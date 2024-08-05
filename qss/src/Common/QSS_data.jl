
 #= """ hold helper datastructures needed for simulation, can be seen as the model in the qss architecture (model-integrator-quantizer)""" =#
struct CommonQSS_data{Z}
    quantum :: Vector{Float64} 
    x :: Vector{Taylor0}  #MVector cannot hold non-isbits
    q :: Vector{Taylor0}
    tx ::  Vector{Float64} 
    tq :: Vector{Float64} 
    d::Vector{Float64} 
    nextStateTime :: Vector{Float64}    
    nextInputTime :: Vector{Float64} # 
    nextEventTime :: MVector{Z,Float64}  
    t::Taylor0# mute taylor var to be used with math functions to represent time
    integratorCache::Taylor0
    taylorOpsCache::Vector{Taylor0}
    finalTime:: Float64   
    savetimeincrement::Float64 
    initialTime :: Float64    
    dQmin ::Float64    
    dQrel ::Float64  
    maxErr ::Float64  
    maxiters ::Int
    savedTimes :: Vector{Vector{Float64}}
    savedVars:: Vector{Vector{Float64}}
    
end
#= ```
data needed only for implicit case
``` =#
struct AexprLiQSS_data{O,M} <: LiQSS_data{O,M}
    cycleDetection::Val{M}
    cacheA::MVector{1,Float64}
    qaux::Vector{MVector{O,Float64}}
    dxaux::Vector{MVector{O,Float64}}
end
struct AmanualLiQSS_data{O,M}<: LiQSS_data{O,M}
    cycleDetection::Val{M}
    a::Vector{Vector{Float64}}
    qaux::Vector{MVector{O,Float64}}
    dxaux::Vector{MVector{O,Float64}}
    olddx::Vector{MVector{O,Float64}}
end
  
struct Options{SU,DU} #single update, double update
    singleUpdate::Val{SU}
    SimulUpdate::Val{DU}
    multiplier::Float64
end

option(su,du,multiplier)=Options(Val(su),Val(du),multiplier)

