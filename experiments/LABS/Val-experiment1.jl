using BenchmarkTools

#=
qss1()=Val(1)
qss2()=Val(2)
qss3()=Val(3)
liqss()=Val(4)

  
function modelSetting(solverType::Val{T}) where {T}
   # println("initializing")
 
    solverType
end
function integrator(solver::Val{T}) where {T}
    computenext(solver)
end
function computenext(::Val{1})
   # println("logic of qss1")
    1
end
function computenext(::Val{2})
    #println("logic of qss2")
    2
end


function test()
    solver=modelSetting(qss1())
    integrator(solver)
end
=#


  
function modelSetting(solverType::Val{T}) where {T}
   # println("initializing")
 
    solverType
end
function integrator(solver::Int) 
    computenext(solver)
end
function computenext(::Val{1})
   # println("logic of qss1")
    1
end
function computenext(::Val{2})
    #println("logic of qss2")
    2
end


function test()
   # solver = Ref(2);
    solver=modelSetting(qss1())
    integrator((solver))
end

#=
qss1()=1
qss2()=2
qss3()=3
liqss()=4

  
function modelSetting(solverType::Int) 
   # println("initializing")
 
    solverType
end
function integrator(solver::Int) 
    computenext(solver)
end
function computenext(solver::Int)
   # println("logic of qss1")
    if solver==1
        1
    elseif solver==2
        2
    else
        3
    end
end



function test()
    solver=modelSetting(qss1())
    integrator(solver)
end
@btime test()   

=#

#=
  
function modelSetting(solverType::String) 
   # println("initializing")
 
    solverType
end
function integrator(solver::String) 
    computenext(solver)
end
function computenext(solver::String)
   # println("logic of qss1")
    if solver=="qss1"
        1
    elseif solver=="qss2"
        2
    else
        3
    end
end



function test()
    solver=modelSetting("qss1")
    integrator(solver)
end

=#
@btime test()  