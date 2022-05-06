struct ModelSettings{N}   
    initConditions ::SVector{N,Float64} 
    jacobian::SMatrix{N,N,Float64}
    finalTime:: Float64   
    initialTime :: Float64    
    dQmin ::Float64    
    dQrel ::Float64  
    solver ::Val{T} where{T}
    savetimeincrement::Float64
end
function ModelSettings(initConditions ,jacobian ,finalTime; initialTime=0.0,dQmin=1e-6,dQrel=1e-3,solver=Val(1),savetimeincrement=0)
    return ModelSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end
#= function ModelSettings(initConditions ,jacobian ,finalTime,solver; initialTime=0.0,dQmin=1e-6,dQrel=1e-3,order=1,savetime=0)
    return ModelSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,order,solver,savetime)
end =#
#= function ModelSettings(initConditions ,jacobian ,finalTime,savetimeincrement; initialTime=0.0,dQmin=1e-6,dQrel=1e-3,solver=Val(1))
    return ModelSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end =#
function ModelSettings(initConditions ,jacobian ,finalTime,savetimeincrement,solver; initialTime=0.0,dQmin=1e-6,dQrel=1e-3)
    return ModelSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end
#saveat(savetimeincrement::Float64)=savetimeincrement
function saveat(savetimeincrement::Float64)
    savetimeincrement
end 
qss1()=Val(1)
qss2()=Val(2)
qss3()=Val(3)
liqss1()=Val(4)
liqss2()=Val(5)










    #=function ModelSettings( it:: Float64,ft:: Float64, dQmin ::Float64, dQrel ::Float64,order :: Int,solver :: Val{T} , initConditions ::Array{Float64} ,jacobian :: Array{Float64, 2} )where {T}
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

