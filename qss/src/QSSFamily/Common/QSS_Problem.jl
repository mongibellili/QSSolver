abstract type qssProblem end
abstract type easyqssProblem end
struct Problem{T,U,D,Z} <: qssProblem
    initConditions ::SVector{T,Float64} 
    jacobian::SVector{T,Function} # if time-variant jac (more efficient to use MMatrix), i suggest you code another modelSetting and specialize when the user calls
    inputVars::SVector{U,Float64} 
    discreteVars::MVector{D,Float64} 
    zcFunctions::SVector{Z,Function} 
    eventHandlerFunctions::SVector{Z,Function} # to each zcFunc there is an eventHandler
    SD::SMatrix{T,T,Int}  # T states can influence  T derivatives
    SZ::SMatrix{T,Z,Int} # T states can influence Z zcfunctions
    HD::SMatrix{D,T,Int}  # D discrete Vars can influence T derivatives
    HZ::SMatrix{D,Z,Int}  # D discrete Vars can influence Z zcfunctions
end
struct TimedProblem{T,U,D,Z}  <: qssProblem 
    initConditions ::SVector{T,Float64} 
    jacobian::SVector{T,Function} 
    inputVars::SVector{U,Function} 
    discreteVars::MVector{D,Float64} 
    zcFunctions::SVector{Z,Function} 
    eventHandlerFunctions::SVector{Z,Function}
    SD::SMatrix{T,T,Int} 
    SZ::SMatrix{T,Z,Int}
    HD::SMatrix{D,T,Int}
    HZ::SMatrix{D,Z,Int}
end
struct EasyProblem{T,U} <: easyqssProblem
    initConditions ::SVector{T,Float64} 
    jacobian::SMatrix{T,T,Float64}  # if time-variant jac (more efficient to use MMatrix), i suggest you code another modelSetting and specialize when the user calls
    inputVars::SVector{U,Float64}  
end
struct EasyTimedProblem{T,U} <: easyqssProblem
    initConditions ::SVector{T,Float64} 
    jacobian::SMatrix{T,T,Float64}  # if time-variant jac (more efficient to use MMatrix), i suggest you code another modelSetting and specialize when the user calls
    inputVars::SVector{U,Function}  
end
function QSS_Problem(initConditions ,jacobian ,inputVars,discreteVars,zcFunctions,eventHandlerFunctions,SD,SZ,HD,HZ)
    try
        return Problem(initConditions, jacobian, inputVars,discreteVars,zcFunctions,eventHandlerFunctions,SD,SZ,HD,HZ)
    catch e
        return TimedProblem(initConditions, jacobian, inputVars,discreteVars,zcFunctions,eventHandlerFunctions,SD,SZ,HD,HZ)
    end
    #if typeof(inputVars)==
    #return ModelSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement,inputVars)#
end
function QSS_Problem(initConditions ,jacobian ,inputVars)    
    try
        return EasyProblem(initConditions, jacobian, inputVars)
    catch
        return EasyTimedProblem(initConditions, jacobian, inputVars)
    end 
end





