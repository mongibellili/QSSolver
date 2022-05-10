
 function updateScheduler(counter,nextStateTime::MVector{T,Float64},nextEventTime :: MVector{Z,Float64} )where{T,Z}   
    
    # which is faster? finding the minimum or implementing a priority queue
    minStateTime=Inf
    minState_index=0  # what if all nextstateTime= Inf ...especially at begining????? min_index stays 0!!!
    minEventTime=Inf
    minEvent_index=0
    ST_STATE=1
    ST_EVENT=2
    for i=1:T
        if nextStateTime[i]<minStateTime
            minStateTime=nextStateTime[i]
            minState_index=i
        end
    end
    for i=1:Z
        if nextEventTime[i] < minEventTime
            minEventTime=nextEventTime[i]
            minEvent_index=i
        end
    end



    if minState_index==0 
        println(" static system! all derivatives are null!")
        return (1,minTime) # later throw exception or maybe draw horizontal lines at initial conditions
    end 
    if minEventTime<minStateTime
       # println("an event N",minEvent_index, "about to occur! at time= ",minEventTime)
        return (minEvent_index,minEventTime,ST_EVENT)
    else
        return (minState_index,minStateTime,ST_STATE)
    end


    
end