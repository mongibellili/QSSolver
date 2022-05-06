#=
#without saving
mutable struct Quantizer    
    computeNext ::Function
    reComputeNext ::Function
    updateQVar ::Function
    solver :: Int    
    function Quantizer(qssSimulator :: QSS_simulator)
        p=new()
        solver=qssSimulator.settings.solver
        if solver =="qss1"
            QSS1_init(p,qssSimulator)
        elseif solver =="qss2"
            #QSS2_init(qssSimulator)
        else 
            #LiQSS_init(qssSimulator)
        end
        p
    end
end
function computeNextTime(solver::Int,quantizer ::Quantizer,index::Int,currentTime::Float64,nextTime::Vector{Float64},x:: Vector{Float64} ,quantum::Vector{Float64})
  #  quantizer.computeNext(quantizer,index,currentTime,nextTime,x,quantum)
  if solver =="qss1"
    QSS1_ComputeNext(quantizer,index,currentTime,nextTime,x,quantum)
  end
end
function reComputeNextTime(solver::Int,quantizer ::Quantizer,index::Int,currentTime::Float64,nextTime::Vector{Float64},x:: Vector{Float64} ,q::Vector{Float64},quantum::Vector{Float64},casesVact::Vector{Float64})
   # quantizer.reComputeNext(quantizer,index,currentTime,nextTime,x,q,quantum,casesVact)
    if solver =="qss1"
   QSS1_reComputeNext(quantizer,index,currentTime,nextTime,x,q,quantum,casesVact)
    end 
end
function updateQ(solver::Int,quantizer ::Quantizer,index::Int,x::Vector{Float64} ,q::Vector{Float64},quantum::Vector{Float64})
    #quantizer.updateQVar(quantizer,index,x,q,quantum)
    if solver =="qss1"
    QSS1_update(quantizer,index,x,q,quantum)
    end
end

=#
#=
mutable struct Quantizer    
    computeNext ::Function
    reComputeNext ::Function
    updateQVar ::Function
    solver :: Int    
    function Quantizer(qssSimulator :: QSS_simulator)
        p=new()
        solver=qssSimulator.settings.solver
        if solver ==1
            QSS1_init(p,qssSimulator)
        elseif solver ==2
            #QSS2_init(qssSimulator)
        else 
            #LiQSS_init(qssSimulator)
        end
        p
    end
end
function computeNextTime(solver::Int,quantizer ::Quantizer,index::Int,currentTime::Float64,nextTime::Vector{Float64},x:: Vector{Array{Float64}} ,quantum::Vector{Float64})
   #quantizer.computeNext(quantizer,index,currentTime,nextTime,x,quantum)
  if solver ==1
    QSS1_ComputeNext(index,currentTime,nextTime,x,quantum)
  end
end
function reComputeNextTime(solver::Int,quantizer ::Quantizer,index::Int,currentTime::Float64,nextTime::Vector{Float64},x:: Vector{Array{Float64}} ,q::Vector{Float64},quantum::Vector{Float64},casesVact::Vector{Float64})
  #  quantizer.reComputeNext(quantizer,index,currentTime,nextTime,x,q,quantum,casesVact)
   if solver ==1
   QSS1_reComputeNext(index,currentTime,nextTime,x,q,quantum,casesVact)
    end 
end
function updateQ(solver::Int,quantizer ::Quantizer,index::Int,x::Vector{Array{Float64}} ,q::Vector{Float64},quantum::Vector{Float64})
   # quantizer.updateQVar(quantizer,index,x,q,quantum)
    if solver ==1
    QSS1_update(index,x,q,quantum)
    end
end
=#