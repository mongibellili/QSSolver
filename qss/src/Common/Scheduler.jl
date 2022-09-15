
 function updateScheduler(nextStateTime::MVector{T,Float64},nextEventTime :: MVector{Z,Float64},nextInputTime :: Vector{Float64} )where{T,Z}   
    
    # which is faster? finding the minimum or implementing a priority queue
    minStateTime=Inf
    minState_index=0  # what if all nextstateTime= Inf ...especially at begining????? min_index stays 0!!!
    minEventTime=Inf
    minEvent_index=0
    minInputTime=Inf
    minInput_index=0
   #=  ST_STATE=1
    ST_INPUT=2
    ST_EVENT=3 =#
    returnedVar=() #  used to print something if something is bad
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
    for i=1:T
        if nextInputTime[i] < minInputTime
            minInputTime=nextInputTime[i]
            minInput_index=i
        end
    end



#=     if minState_index==0 
        println(" static system! all derivatives are null!")
        return (1,minTime) # later throw exception or maybe draw horizontal lines at initial conditions
    end  =#
    if minEventTime<minStateTime
       # println("an event N",minEvent_index, "about to occur! at time= ",minEventTime)
       if minInputTime<minEventTime
          
         #return (minInput_index,minInputTime,:ST_INPUT)
         returnedVar=(minInput_index,minInputTime,:ST_INPUT)
       else
       
       
       #return (minEvent_index,minEventTime,:ST_EVENT)
       returnedVar=(minEvent_index,minEventTime,:ST_EVENT)
       
       end
    else
       
        if minInputTime<minStateTime
            # println("an event N",minEvent_index, "about to occur! at time= ",minEventTime)
            # return (minInput_index,minInputTime,:ST_INPUT)
            returnedVar=(minInput_index,minInputTime,:ST_INPUT)
         else
       
             #return (minState_index,minStateTime,:ST_STATE)
             returnedVar=(minState_index,minStateTime,:ST_STATE)
         end
    end
    if returnedVar[1]==0
        println("scheduler *******nextEventTime= ",nextEventTime)
        println("scheduler *******nextInputTime= ",nextInputTime)
        println("scheduler *******nextstateTime= ",nextStateTime)
        println("also the whole system may be static! developer: fill remaing points with same values and exit!")
    end
    return returnedVar


    
end

#= 
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

function updateScheduler(nextStateTime::MVector{T,Float64},nextInputTime :: Vector{Float64} )where{T}   
      
    # which is faster? finding the minimum or implementing a priority queue
    minStateTime=Inf
    minState_index=0  # what if all nextstateTime= Inf ...especially at begining????? min_index stays 0!!!
    minInputTime=Inf
    minInput_index=0
    ST_STATE=1
    ST_INPUT=2
    for i=1:T
        if nextStateTime[i]<minStateTime
            minStateTime=nextStateTime[i]
            minState_index=i
        end
    end
    for i=1:T
        if nextInputTime[i] < minInputTime
            minInputTime=nextInputTime[i]
            minInput_index=i
        end
    end



#=     if minState_index==0 
        println(" static system! all derivatives are null!")
        return (1,minTime) # later throw exception or maybe draw horizontal lines at initial conditions
    end  =#
    if minInputTime<minStateTime
       # println("an event N",minEvent_index, "about to occur! at time= ",minEventTime)
        return (minInput_index,minInputTime,:ST_INPUT)
    else
        return (minState_index,minStateTime,:ST_STATE)
    end


  end   =#