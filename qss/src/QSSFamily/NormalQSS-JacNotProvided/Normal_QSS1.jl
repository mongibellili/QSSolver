using InteractiveUtils

function computeDerivative(index::Int, ::Val{1}   ,x::MVector{O,Float64} , q::MVector{T,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64}  ) where{T,O}   
  
   x[(2)*index]=Main.f(Val(index),q)
 # x[(2)*index]=funcX(Val(index),q,tq,1)
#=    if index==1
    x[(2)*index]=q[2]
  else
    x[(2)*index]=-q[1]-q[2]+1
  end  =#

#=   if index==1
    x[(2)*index]=Main.funcX1(q,tq,1)
  else
    x[(2)*index]=Main.funcX2(q,tq,1)
  end =#

#=   f = Symbol(:funcX,index)
  x[(2)*index]=@eval $f($q,$tq,1)  =#
  #x[(2)*index]=getfield(Main, Symbol(:funcX,index))(q,tq,1)
  #s= Symbol(:funcX,index)
  #f = getfield(Main, Symbol(:funcX,index))
  #display(@code_lowered x[(2)*index]=getfield(Main, Symbol(:funcX,index))(q,tq,1))
  #display(@code_warntype x[(2)*index]=getfield(Main, Symbol(:funcX,index))(q,tq,1))
  #f=nothing
 #display(@code_lowered jacobian[1](q,tq,1) )
  # display(@code_typed jacobian[1](q,tq,1)  )
   #display(@code_llvm jacobian[1](q,tq,1))
  # println("-------------------------------------")
   #display(@code_llvm tempF(q,tq,1) )
   #display(@code_lowered tempF(q,tq,1) )
  # display(@code_typed tempF(q,tq,1) )
  #---------@code_warntype-----------
  
 end