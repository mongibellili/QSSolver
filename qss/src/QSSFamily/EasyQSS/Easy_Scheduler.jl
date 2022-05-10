
 function updateScheduler(nextStateTime::MVector{T,Float64} )where{T}   
    minTime=Inf
    min_index=0  # what if all nextstateTime= Inf ...especially at begining????? min_index stays 0!!!

    for i=1:T
        if nextStateTime[i]<minTime
            minTime=nextStateTime[i]
            min_index=i
        end
    end
#=     mintimevalue[1]=minTime
    minIndex[1]=min_index =#

    if min_index==0
        println(" static system! all derivatives are null!")
        return (1,minTime) # later throw exception or maybe draw horizontal lines at initial conditions
    end 
    (min_index,minTime)
end