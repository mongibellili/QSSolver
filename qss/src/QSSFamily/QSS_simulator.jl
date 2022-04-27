# struct that glues the other 3 structs settings, data, time, model
# information flow is done through it
struct QSS_simulator
   
    x :: Vector{Array{Float64}}
    tx ::  Vector{Array{Float64}} 
   #= q :: MVector{L,Float64} where {L} 
    tq :: MVector{M,Float64} where {M}
    quantum :: MVector{N,Float64} where {N} 
    nextStateTime :: MVector{P,Float64} where {P}  
    =#
    q :: Vector{Float64} 
    tq :: Vector{Float64} 
    quantum :: Vector{Float64} 
    nextStateTime :: Vector{Float64} 
    wRegister::MVector{2,Float64}(undef)#time,minValue   
    minIndex :: MVector{1,Float64}(undef) 
    jacobian :: Array{Float64, 2}   
    dep ::  Vector{Array{Int}}  



end
function QSS_simulate(settings :: ModelSettings)
    #println("........started simulation.........")

        const d.states=states
        const d.order=order
        #create Vectors
        d.quantum = Vector{Float64}(undef, states)
        d.x =  Vector{Array{Float64}}(undef, (order+1)*states)#case1 for x case2 for DERx...
       # d.q =  Vector{Float64}(undef, (order+1)*states)
       d.q = Vector{Float64}(undef, (order+1)*states)
        #create arrays inside Vectors
        for i = 1:(order+1)*states 
            d.x[i]=Array{Float64}[]
           # d.q[i]=Array{Float64}[]
        end
     
        d
    
   QSS_integrate(qssSimulator);
end
