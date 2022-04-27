
function computeDerivative(states, index, order, jacobian , x, q, tx, tq)
  #for now i will not consider time ie derx=f(x)   
  x[(order+1)*index]=0
  for j = 1:states
    x[(order+1)*index] += jacobian[j,index] * q[(order+1)*j-order] # use jacobian[j,index] when jacobian=transpose(jacobian)
    #x[(order+1)*index] += jacobian[index,j] * q[(order+1)*j-order] 
  end
end

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

#=function computeDerivative(states::Int, index::Int, order::Int, jacobian::SMatrix{2,2,Float64} , x::MVector{4,Float64}, q::MVector{4,Float64}, tx::MVector{2,Float64}, tq::MVector{2,Float64})
  #for now i will not consider time ie derx=f(x)   
  #der = 0
x[(order+1)*index]=0
  for j = 1:states
    x[(order+1)*index] += jacobian[index, j] * q[(order+1)*j-order]
  end
  #x[(order+1)*index]= der


end=#
