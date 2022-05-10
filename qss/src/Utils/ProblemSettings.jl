struct ProblemSetting
    finalTime:: Float64   
    initialTime :: Float64    
    dQmin ::Float64    
    dQrel ::Float64  
    solver ::Val{V} where{V}
    savetimeincrement::Float64 
end
# the user might enter these in a different order: no compile error, but semantic error
function ProblemSettings(finalTime,savetimeincrement,solver, initialTime,dQmin,dQrel) 
    return ProblemSetting(finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end
function ProblemSettings(finalTime,savetimeincrement,solver,dQmin,dQrel; initialTime=0.0) 
    return ProblemSetting(finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end
function ProblemSettings(finalTime,savetimeincrement,solver; initialTime=0.0,dQmin=1e-6,dQrel=1e-3) 
    return ProblemSetting(finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end
function ProblemSettings(finalTime,savetimeincrement;solver=Val(1), initialTime=0.0,dQmin=1e-6,dQrel=1e-3) 
    return ProblemSetting(finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end
function ProblemSettings(finalTime;savetimeincrement=0.1,solver=Val(1), initialTime=0.0,dQmin=1e-6,dQrel=1e-3) 
    return ProblemSetting(finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end

function saveat(savetimeincrement::Float64)
    savetimeincrement
end 
qss1()=Val(1)
qss2()=Val(2)
qss3()=Val(3)
liqss1()=Val(4)
liqss2()=Val(5)










    #=function ProblemSettings( it:: Float64,ft:: Float64, dQmin ::Float64, dQrel ::Float64,order :: Int,solver :: Val{T} , initConditions ::Array{Float64} ,jacobian :: Array{Float64, 2} )where {T}
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

#= function ProblemSettings(initConditions ,jacobian ,finalTime,solver; initialTime=0.0,dQmin=1e-6,dQrel=1e-3,order=1,savetime=0)
    return ProblemSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,order,solver,savetime)
end =#
#= function ProblemSettings(initConditions ,jacobian ,finalTime,savetimeincrement; initialTime=0.0,dQmin=1e-6,dQrel=1e-3,solver=Val(1))
    return ProblemSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement)
end =#