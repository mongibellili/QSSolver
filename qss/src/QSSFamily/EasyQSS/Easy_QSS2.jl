function integrateState(j::Int,::Val{2},elapsed::Float64,x::MVector{O,Float64})where{O} 
  x[3*j-2]=x[3*j-2]+ elapsed * x[3*j-1] +0.5*elapsed*elapsed*x[3*j] 
  x[3*j-1]=x[3*j-1]+ elapsed * x[3*j] 
  #println("xi",j,"= ",x[3*j-2])
  #println("derxi",j,"= ",x[3*j-1])
end
function reIntegrateState(j::Int,::Val{2},elapsed::Float64,elapsedQ::Float64,x::MVector{O,Float64}, q::MVector{R,Float64})where{O,R} 
  x[3*j-2]=x[3*j-2]+ elapsed * x[3*j-1] +0.5*elapsed*elapsed*x[3*j] 
  #where should I update q=q+e*derq????
  q[2*j-1]= x[3*j-2] # 
  #println("x",j,"= ",x[3*j-2])
end




function computeInitDerivative(j::Int, ::Val{2}, jacobian::SVector{T,SVector{T,Float64}} ,x::MVector{O,Float64} , q::MVector{R,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64} , inputVars::SVector{U,Float64}  )where{R,T,O,U}
  #for now i will not consider time ie derx=f(x)   
  x[(3)*j-1]=0
  for k = 1:T
    x[(3)*j-1] += jacobian[j][k] * q[2*k-1]      
  end
  x[(3)*j-1] +=  inputVars[j] 
end
  function computeInitsecondDerivative(j::Int, ::Val{2}, jacobian::SVector{T,SVector{T,Float64}} ,x::MVector{O,Float64} , q::MVector{R,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64}, inputVars::SVector{U,Float64}  )where{R,T,O,U}
    #for now i will not consider time ie derx=f(x)    
    x[3*j] =0
    for k = 1:T
     #x[(3)*j] += derjacobian[j][k] * q[2*k-1]  
   
     x[3*j] +=jacobian[j][k] * q[2*k] 
    end
end


function computeDerivative(t::Float64,index::Int,j::Int, ::Val{2}, jacobian::SVector{T,SVector{T,Float64}}  ,x::MVector{O,Float64} , q::MVector{R,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64} ,inputVars::SVector{U,Float64}  )where{R,T,O,U}
  #for now i will not consider time ie derx=f(x)   
  x[(3)*j-1]=0
  x[3*j] =0
  for k = 1:T               
    x[(3)*j-1] += jacobian[j][k] * q[2*k-1]  # this q should be updated using q=q+elapsedQ*derQ if not already been updated using q=x ;(qindex, and qj those that undergone reintegrateState)
    x[3*j] +=jacobian[j][k] *  q[(2)*k]
    
  end 
  x[(3)*j-1] +=  inputVars[j] 
 # println("derx",j,"= ",x[(3)*j-1])
end

function updateQ(::Val{2},index::Int, x::MVector{O,Float64}, q::MVector{R,Float64}, quantum::MVector{T,Float64}) where{R,T,O}
  q[(2)*index-1] = (x[(3)*index-2]) #q=x 
  q[(2)*index] = (x[(3)*index-1]) 
end
function updatederQ(::Val{2},index::Int, x::MVector{O,Float64}, q::MVector{R,Float64}, quantum::MVector{T,Float64}) where{R,T,O}
  q[(2)*index] = (x[(3)*index-1])  #derq= derx
  #println("updatedDerQ",index,"= ",q[2*index])
end
function computeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::MVector{O,Float64}, quantum::MVector{T,Float64})where{T,O}
  if (x[(3)*index]) != 0
    nextTime[index] = currentTime +sqrt( abs(quantum[index] / (x[(3)*index])))
  else
    nextTime[index] = Inf
  end
end
function reComputeNextTime(::Val{2}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::MVector{O,Float64}, q::MVector{R,Float64}, quantum::MVector{T,Float64})where{T,O,R}
  coef=@SVector [q[2*index-1] - x[3*index-2] - quantum[index],q[2*index] - x[3*index-1] , -x[3*index]]#[c,b,a] coef of quadratic eq
  time1 = currentTime + minPosRoot(coef, Val(2))
  coef=setindex(coef,q[2*index-1] - x[3*index-2] + quantum[index],1)
  time2 = currentTime + minPosRoot(coef, Val(2))
  nextTime[index] = time1 < time2 ? time1 : time2 
end
 
