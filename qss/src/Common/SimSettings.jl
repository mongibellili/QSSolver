#struct to hold simulation run setting...solver setting










    #=function SimSettings( it:: Float64,ft:: Float64, dQmin ::Float64, dQrel ::Float64,order :: Int,solver :: Val{T} , initConditions ::Array{Float64} ,jacobian :: Array{Float64, 2} )where {T}
        mSettings = new()
        mSettings.initConditions = initConditions
        mSettings.jacobian = jacobian
        mSettings.initialTime = it
        mSettings.finalTime = ft*1.0
        mSettings.dQmin = dQmin
        mSettings.dQrel = dQrel
        mSettings.order=order
        mSettings.solver = solver
        mSettings.states= computeStates(jacobian)
    
        mSettings
    end
    
    
    function computeStates(jacobian :: Array{Float64, 2})
        # return number of rows
        return size(jacobian,1)
    end
    =#

#= function SimSettings(initConditions ,jacobian ,finalTime,solver; initialTime=0.0,dQmin=1e-6,dQrel=1e-3,order=1,savetime=0)
    return SimSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,order,solver,savetime)
end =#
#= function SimSettings(initConditions ,jacobian ,finalTime,savetimeincrement; initialTime=0.0,dQmin=1e-6,dQrel=1e-3,solver=Val(1))
    return SimSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end =#