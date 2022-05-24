#= function integrateState(j::Int,::Val{1},elapsed::Float64,x::MVector{O,Float64})where{O} 
  x[2*j-1]=x[2*j-1]+ elapsed * x[2*j] 
end =#
#= function reIntegrateState(j::Int,::Val{1},elapsed::Float64,x::MVector{O,Float64})where{O} 
  x[2*j-1]=x[2*j-1]+ elapsed * x[2*j] 
end =#
function computeDerivative(index::Int, ::Val{1} , jacobian::SVector{T,Function}  ,x::MVector{O,Float64} , q::MVector{T,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64},  d::MVector{D,Float64}  )where{T,O,D}   
    x[(2)*index] = jacobian[index](q,d,tq,1)        
end



function computeNextEventTime(elapsed,counter,solver,j,zcFunctions,oldsignValue,dt,currentTime,  nextEventTime, q,d,tq,x, quantum) #later specify args
  #logic of next event time
  if 0.05555>currentTime>0.0555
    counter[1]=5
  end

  if 0.9999>currentTime>0.999
    counter[1]=5
  end
  output=zcFunctions[j](q,d,tq,1) # 1 for order 1

 # println("output2= ",output[2])
  if oldsignValue[j,1] != output[1]
    nextEventTime[j]=currentTime 
    if counter[1]>0
      println("old sign of",j,"  =  ",oldsignValue[j,1])
      println("new sign of",j,"  =  ",output[1])
      println("q  =  ",q)
      println("x  =  ",x)
      println("current time= ",currentTime)
      println("schedule event now at= ",nextEventTime[j])
    end
    counter[1]-=1
  else
    nextEventTime[j] = Inf  #later we can estimate the time
  
    
  #  println("newoutput2= ",newOutput[2])
    coef2=-(output[2]-oldsignValue[j,2])/elapsed  # der(expr)
    coef1=oldsignValue[j,2] # oldExpr value; newExpr= oldExpr+ desiredTime * der(expr)=0
    #= coef=@SVector [coef1, coef2]
    time1 = currentTime + minPosRoot(coef, Val(1))
    q[1]=q[1]-2*dt # try q=q-dt
    newOutput=zcFunctions[j](q,d,tq[j]+dt,1)
   # println("newoutput22= ",newOutput[2])
    coef2=-(newOutput[2]-output[2])/dt  # der(expr)
    coef=setindex(coef,coef2,2)
    time2 = currentTime + minPosRoot(coef, Val(1))
    # i think if both time 1 and time 2 =infinite then derx=0...I doubt the mpr be negative for both cases
    nextEventTime[j] = time1 < time2 ? time1 : time2  =#
#=     if counter[1]>0
     println("schedule event",j," at= ",nextEventTime[j])
    end =#
#=     if (coef2) != 0
      nextEventTime[j] = currentTime +abs(coef1 /(coef2))
      #println("schedule event at= ",nextEventTime[j])
    else
      nextEventTime[j] = Inf
      #println("schedule event at Inf")
    end   =#  
  end
  oldsignValue[j,1]=output[1]
  oldsignValue[j,2]=output[2]
  counter[1]-=1
end
#= function computeDerivative(index::Int, ::Val{1} , jacobian::SVector{T,Function} ,x::MVector{O,Float64} , q::MVector{T,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64},  inputVars::SVector{U,Function}  )where{T,O,U}
  #function computeDerivative(index::Int, order::Int, jacobian::SMatrix{T,T,Float64} ,x::MVector{O,Float64} , q::MVector{O,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64})where{T,O}
    #for now i will not consider time ie derx=f(x)   
    x[(2)*index]=0
    for j = 1:T
     # x[(order+1)*index] += jacobian[j,index] * q[(order+1)*j-order] # use jacobian[j,index] when jacobian=transpose(jacobian)
     # x[(order+1)*index] += jacobian[index+(j-1)*T] * q[(order+1)*j-order] 
      x[(2)*index] += jacobian[index][j] * q[j]   # later try x[]=muladd(jacobian[index][j] , q[j],x[])
    end
       x[(2)*index] +=  inputVars[index](tx[index])   # later try x[]=muladd(jacobian[index][j] , q[j],x[]) # when number of u is less than 
  end =#
#= function updateQ(::Val{1},index::Int, x::MVector{O,Float64}, q::MVector{T,Float64}, quantum::MVector{T,Float64}) where{T,O}  
  q[index] = (x[(2)*index-1])#q=x
end
function computeNextTime(::Val{1}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::MVector{O,Float64}, quantum::MVector{T,Float64})where{T,O}
  if (x[(2)*index]) != 0
    nextTime[index] = currentTime + abs(quantum[index] / (x[(2)*index]))
  else
    nextTime[index] = Inf
  end
end
function reComputeNextTime(::Val{1}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::MVector{O,Float64}, q::MVector{T,Float64}, quantum::MVector{T,Float64})where{T,O}
  coef=@SVector [q[index] - (x[(2)*index-1]) - quantum[index], -x[(2)*index]]
  time1 = currentTime + minPosRoot(coef, Val(1))
  coef=setindex(coef,q[index] - (x[(2)*index-1]) + quantum[index],1)
  time2 = currentTime + minPosRoot(coef, Val(1))
  # i think if both time 1 and time 2 =infinite then derx=0...I doubt the mpr be negative for both cases
  nextTime[index] = time1 < time2 ? time1 : time2 
  
  end
  =#

  #= mpr=-1
if x[(2)*index] == 0 
  mpr = Inf
else 
    mpr1 = (q[(2)*index-1] - (x[(2)*index-1]) - quantum[index]) / x[(2)*index];
    mpr2 = (q[(2)*index-1] - (x[(2)*index-1]) + quantum[index]) / x[(2)*index];            
    if mpr1 < 0
        mpr1 = Inf
    end
    if mpr2 < 0
      mpr2 = Inf
    end
    mpr=mpr1 < mpr2 ? mpr1 : mpr2
end
nextTime[index] = currentTime +mpr
end =#

#mpr=-1
# epselon=10
#=if abs(x[(2)*index])< epselon 
  mpr = Inf =#
 # println("inside recompute qss1")
#=  if x[(2)*index] == 0 
  mpr = Inf
else 
    mpr1 = (q[index] - (x[(2)*index-1]) - quantum[index]) / x[(2)*index];
    mpr2 = (q[index] - (x[(2)*index-1]) + quantum[index]) / x[(2)*index];            
    if mpr1 < 0
        mpr1 = Inf    
      if mpr2 < 0
        mpr = Inf
      else
        mpr=mpr2
      end
    else
      if mpr2 < 0
        mpr = mpr1
      else
        mpr=mpr1 < mpr2 ? mpr1 : mpr2
      end
    end
end
nextTime[index] = currentTime +mpr =#






 




