
function computeDerivative(index::Int, ::Val{1}   ,x::MVector{O,Float64} , q::MVector{T,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64}, jacobian::SVector{T,Function}   ) where{T,O}   
   x[(2)*index] = jacobian[index](q,tq,1)       # allocates like hell 
   #x[(2)*index]=Main.f(Val(index),q)
  # x[(2)*index]=funcX(Val(index),q,tq,1)
#=    if index==1
    x[(2)*index]=q[2]
  else
    x[(2)*index]=-q[1]-q[2]
  end  =#
 end