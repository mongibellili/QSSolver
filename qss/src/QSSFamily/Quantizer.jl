mutable struct Quantizer    
    computeNext ::Function
    reComputeNext ::Function
    updateQVar ::Function
    solver :: String    
    function Quantizer(qssSimulator :: QSS_simulator)
        p=new()
        solver=qssSimulator.settings.solver
        if solver =="qss1"
            QSS1_init(p,qssSimulator)
        elseif solver =="qss2"
            #QSS2_init(qssSimulator)
        else 
            #Liqss_init(qssSimulator)
        end
        p
    end
end

function computeNextTime(quantizer ::Quantizer,index::Int,currentTime::Float64,nextTime::Vector{Float64},x::Vector{Array{Float64}} ,quantum::Vector{Float64})
    quantizer.computeNext(quantizer,index,currentTime,nextTime,x,quantum)
end

function reComputeNextTime(quantizer ::Quantizer,index::Int,currentTime::Float64,nextTime::Vector{Float64},x::Vector{Array{Float64}} ,q::Vector{Array{Float64}},quantum::Vector{Float64})
    quantizer.reComputeNext(quantizer,index,currentTime,nextTime,x,q,quantum)
end





function updateQ(quantizer ::Quantizer,index::Int,x::Vector{Array{Float64}} ,q::Vector{Array{Float64}},quantum::Vector{Float64})
    quantizer.updateQVar(quantizer,index,x,q,quantum)
end
