 #=
 #without saving
 function computeDerivative(states ::Int,index ::Int,order::Int,qssmodel :: QSSmodel,x ::  Vector{Float64} ,q::Vector{Float64} ,tx:: Vector{Float64} ,tq::Vector{Float64} )
    #for now i will not consider time ie derx=f(x)   
     der=0  
   # if length(x[(order+1)*index])!=0
    #       pop!(x[(order+1)*index])  
    #end
    for j = 1:states
       
        #lastQ=last(q[(order+1)*j-order])
        #der+=matrixEntry*lastQ
          der+=qssmodel.jacobian[index,j] * q[(order+1)*j-order]
    end
      x[(order+1)*index]=der   

    
end
=#

function computeDerivative(states ::Int,index ::Int,order::Int,qssmodel :: QSSmodel,x ::  Vector{Array{Float64}} ,q::Vector{Float64} ,tx:: Vector{Array{Float64}} ,tq::Vector{Float64} )
  #for now i will not consider time ie derx=f(x)   
   der=0  

  for j = 1:states
     
      #lastQ=last(q[(order+1)*j-order])
      #der+=matrixEntry*lastQ
        der+=qssmodel.jacobian[index,j] * q[(order+1)*j-order]
  end
    push!(x[(order+1)*index],der)   

  
end