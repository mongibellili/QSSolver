function QSS1_init(quantizer ::Quantizer,qssSimulator ::QSS_simulator)
    quantizer.computeNext=QSS1_ComputeNext
    quantizer.updateQVar=QSS1_update
end
function QSS1_ComputeNext(quantizer::Quantizer,index::Int,currentTime::Float64,nextTime::Vector{Float64},x::Vector{Array{Float64}} ,quantum::Vector{Float64})
    println("qss1_computeNext")
    if last(x[(2)*index])!=0
        nextTime[index]=currentTime+abs(quantum[index]/last(x[(2)*index]))
    else
        nextTime[index]=Inf
    end
end
function QSS1_update(quantizer ::Quantizer,index::Int,x::Vector{Array{Float64}} ,q::Vector{Array{Float64}},quantum::Vector{Float64})
    #q=x
    push!(q[(2)*index],last(x[(2)*index]))# not tested
end