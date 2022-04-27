
 function updateScheduler(states::Int,nextStateTime::MVector{2,Float64},mintimevalue::MVector{1,Float64},minIndex::MVector{1,Int})
    
    minTime=Inf
    min_index=0
    for i=1:states
        if nextStateTime[i]<minTime
            minTime=nextStateTime[i]
            min_index=i
        end
    end
    mintimevalue[1]=minTime
    minIndex[1]=min_index


end