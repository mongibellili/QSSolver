function computeInitialDerivative(states ::Int,index ::Int,order::Int,qssmodel :: QSSmodel,x :: Vector{Array{Float64}} ,q::Vector{Array{Float64}} ,tx::Vector{Array{Float64}} ,tq::Vector{Array{Float64}} )
    #for now i will not consider time ie derx=f(x)   
    der=0   
    for j = 1:states
       
        #lastQ=last(q[(order+1)*j-order])
        #der+=matrixEntry*lastQ
        der+=qssmodel.jacobian[index,j] * last(q[(order+1)*j-order])
    end
    push!(x[(order+1)*index],der)    
    
end