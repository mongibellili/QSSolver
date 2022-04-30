function updateQ(::Val{1},index::Int, x::MVector{O,Float64}, q::MVector{O,Float64}, quantum::MVector{T,Float64}) where{T,O}
  #q=x
  q[(2)*index-1] = (x[(2)*index-1])
  #q[index]=1.5
end
function computeNextTime(::Val{1}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::MVector{O,Float64}, quantum::MVector{T,Float64})where{T,O}
  if (x[(2)*index]) != 0
    # if 0!=0
    nextTime[index] = currentTime + abs(quantum[index] / (x[(2)*index]))
  else
    nextTime[index] = Inf
  end
end
function reComputeNextTime(::Val{1}, index::Int, currentTime::Float64, nextTime::MVector{T,Float64}, x::MVector{O,Float64}, q::MVector{O,Float64}, quantum::MVector{T,Float64})where{T,O}
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
 if x[(2)*index] == 0 
  mpr = Inf
else 
    mpr1 = (q[(2)*index-1] - (x[(2)*index-1]) - quantum[index]) / x[(2)*index];
    mpr2 = (q[(2)*index-1] - (x[(2)*index-1]) + quantum[index]) / x[(2)*index];            
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
nextTime[index] = currentTime +mpr
end
 




#=   coef=@SVector [q[(2)*index-1] - (x[(2)*index-1]) - quantum[index], -x[(2)*index]]
  time1 = currentTime + minPosRoot(coef, Val(1))# 1 because order1
  coef=setindex(coef,q[(2)*index-1] - (x[(2)*index-1]) + quantum[index],1)
  time2 = currentTime + minPosRoot(coef, Val(1))# 1 because order1
  #println("time1= ",time1);println("time2= ",time2)
  nextTime[index] = time1 < time2 ? time1 : time2 =#