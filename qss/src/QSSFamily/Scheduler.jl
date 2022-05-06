
 function updateScheduler(nextStateTime::MVector{T,Float64})where{T}
    
    minTime=Inf
    min_index=0  # what if all nextstateTime= Inf ????? min_index stays 0!!!
    for i=1:T
        if nextStateTime[i]<minTime
            minTime=nextStateTime[i]
            min_index=i
        end
    end
#=     mintimevalue[1]=minTime
    minIndex[1]=min_index =#

 (min_index,minTime)
end