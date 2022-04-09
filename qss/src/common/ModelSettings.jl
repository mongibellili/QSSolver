mutable struct ModelSettings    
    initConditions ::Array{Float64} 
    jacobian :: Array{Float64, 2}   
    initialTime :: Float64  
    finalTime:: Float64   
    dQmin ::Float64    
    dQrel ::Float64
    order :: Int  
    states :: Int       
    solver :: String 
  

    function ModelSettings( it:: Float64,ft:: Float64, dQmin ::Float64, dQrel ::Float64,order :: Int,solver :: String , initConditions ::Array{Float64} ,jacobian :: Array{Float64, 2} )
        mSettings = new()
        mSettings.initConditions = initConditions
        mSettings.jacobian = jacobian
        mSettings.initialTime = it
        mSettings.finalTime = ft
        mSettings.dQmin = dQmin
        mSettings.dQrel = dQrel
        mSettings.order=order
        mSettings.solver = solver
        mSettings.states= computeStates(jacobian)
    
        mSettings
    end
    function computeStates(jacobian :: Array{Float64, 2})
        # dimension of jac
        return 2
    end
end