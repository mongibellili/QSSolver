mutable struct QSSdata    
    quantum :: Vector{Float64}   
    x :: Vector{Array{Float64}} 
    q :: Vector{Array{Float64}} 
    states :: Int 
    function QSSdata(states :: Int, initConditions ::Array{Float64} )
        d = new()
        d.states=states
        d.quantum = Vector{Float64}(undef, states)
        d.x = Vector{Array{Float64}}(undef, states)
        d.q = Vector{Array{Float64}}(undef, states)
        for i in 1:states
            d.x[i]=Array{Float64}[]
            push!(d.x[i],initConditions[i])#push===fill in first element of arr
            d.q[i]=Array{Float64}[]
            push!(d.q[i],initConditions[i])
        end
        d
    end
end
#given an initialTime and number of states
# it initializes t for all x and all q
mutable struct QSStime  
    time :: Float64
    nextStateTime :: Vector{Float64}   
    tx :: Vector{Array{Float64}} 
    tq :: Vector{Array{Float64}} 
    states :: Int 
    minValue :: Float64   
    minIndex :: Int 
    function QSStime(states :: Int , initTime :: Float64)
        t = new()
        t.states=states
        t.time= initTime
        t.nextStateTime=Vector{Float64}(undef, states)
        t.tx = Vector{Array{Float64}}(undef, states)
        t.tq = Vector{Array{Float64}}(undef, states)
        for i in 1:states
            t.tx[i]=Array{Float64}[]
            push!(t.tx[i],initTime)
            t.tq[i]=Array{Float64}[]
            push!(t.tq[i],initTime)
        end
        t
    end
end
mutable struct QSSmodel 
    
    jacobian :: Array{Float64, 2}   
    dep ::  Array{Float64, 2}   
   
    function QSSmodel(jac :: Array{Float64, 2}  )
        m = new()
        m.jacobian= jac
        m.dep = createDependency(jac)
        m
    end
    function createDependency(jac :: Array{Float64, 2}  )
       # extract dependency matrix from jac
        return jac 
    end

end