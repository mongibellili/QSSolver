function QSS1_init(quantizer ::Quantizer,qssSimulator ::QSS_simulator)
    quantizer.computeNext=QSS1_ComputeNext
    quantizer.reComputeNext=QSS1_reComputeNext
    quantizer.updateQVar=QSS1_update
end
  function QSS1_ComputeNext(quantizer::Quantizer,index::Int,currentTime::Float64,nextTime::Vector{Float64},x:: Vector{Float64} ,quantum::Vector{Float64})
    
    if (x[(2)*index])!=0
          nextTime[index]=currentTime+abs(quantum[index]/(x[(2)*index]))
    else
          nextTime[index]=Inf
    end
end
  function QSS1_reComputeNext(quantizer::Quantizer,index::Int,currentTime::Float64,nextTime::Vector{Float64},x:: Vector{Float64},q::Vector{Float64} ,quantum::Vector{Float64})
    
  
    coef = Vector{Float64}(undef, 2)
      coef[1]=q[(2)*index-1]-(x[(2)*index-1])-quantum[index] #q-x-deltaQ
      coef[2]=-(x[(2)*index])  #-derX
    time1=currentTime+minPosRoot(coef,1)# 1 because order1
      coef[1]=q[(2)*index-1]-(x[(2)*index-1])+ quantum[index] #q-x-deltaQ
    time2=currentTime+minPosRoot(coef,1)# 1 because order1
    #println("time1= ",time1);println("time2= ",time2)
    nextTime[index]= time1 < time2 ? time1 : time2




end
  function QSS1_update(quantizer ::Quantizer,index::Int,x::Vector{Float64} ,q::Vector{Float64},quantum::Vector{Float64})
    #q=x
   #= if length(q[(2)*index-1])!=0
        pop!(q[(2)*index-1])  
    end
    push!(q[(2)*index-1],(x[(2)*index-1]))# not tested
    =#
      q[(2)*index-1]=x[(2)*index-1]
end