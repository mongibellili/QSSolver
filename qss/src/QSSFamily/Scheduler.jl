function updateScheduler(qsstime::QSStime)
    states=qsstime.states
    minTime=Inf
    min_index=0
    for i=1:states
        if qsstime.nextStateTime[i]<minTime
            minTime=qsstime.nextStateTime[i]
            min_index=i
        end
    end
    qsstime.time=minTime
    qsstime.minIndex=min_index

end