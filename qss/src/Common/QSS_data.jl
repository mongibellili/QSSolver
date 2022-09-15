#hold helper datastructures needed for simulation


struct QSS_data{T,Z}
    quantum :: Vector{Float64} 
 #=    x :: MVector{O,Float64}
    q :: MVector{R,Float64}  =#
    x :: Vector{Taylor0{Float64}}
  #  derx :: Vector{Taylor0{Float64}}# used as a cache derx=f...did not increase performance as f returned in cache1 
    q :: Vector{Taylor0{Float64}}
    tx ::  MVector{T,Float64} 
    tq :: MVector{T,Float64} 
    nextStateTime :: MVector{T,Float64}    
    nextInputTime :: Vector{Float64}  
    nextEventTime :: MVector{Z,Float64}  
    t::Taylor0{Float64}
    integratorCache::Taylor0{Float64}
    order::Int
    savedVars :: Vector{Array{Taylor0{Float64}}} #has to be vector (not SA) cuz to be resized in integrator
    savedTimes :: Vector{Float64}  
    taylorOpsCache::Vector{Taylor0{Float64}}

    finalTime:: Float64   
    savetimeincrement::Float64 
    initialTime :: Float64    
    dQmin ::Float64    
    dQrel ::Float64  
end

#= 
function qss_Data1(initConditions::SVector{T,Float64},discreteVars::Vector{Float64},quantum:: Vector{Float64} ,tx::  MVector{T,Float64} ,tq::  MVector{T,Float64} ,nextStateTime::  MVector{T,Float64} ,nextInputTime:: Vector{Float64} ,nextEventTime::  MVector{Z,Float64}  ,cacheSize::Int,SD::SVector{T,SVector{T,Int}},SZ::SVector{T,SVector{Z,Int}},HD::SVector{D,SVector{T,Int}},HZ::SVector{D,SVector{Z,Int}},ZC_jacobian::SVector{Z,SVector{T,Basic}})where{T,Z,D}
    return Static_QSS_data{T,Z,D}(initConditions,discreteVars,quantum,tx,tq,nextStateTime,nextInputTime ,nextEventTime ,cacheSize,SD,SZ,HD,HZ,ZC_jacobian)
end =#  #precise saving goes with this approach

function saveat(savetimeincrement::Float64)
    savetimeincrement
end 
qss1()=Val(1)
qss2()=Val(2)
qss3()=Val(3)
liqss1()=Val(4)
liqss2()=Val(5)


#= getOrderfromSolverMethod(::Val{1})=1
getOrderfromSolverMethod(::Val{2})=2
getOrderfromSolverMethod(::Val{3})=3 =#
#getOrderfromSolverMethod(::Val{4})=1 #liqss1 #having more than 3 methods is causing an error 
#getOrderfromSolverMethod(::Val{5})=2  #liqss2
function getOrderfromSolverMethod(::Val{V}) where {V}  # @generated and inline did not enhance performance
    
    if V==1 || V==2 || V==3
        return V
    else
        return V-3
    end
end


#= struct Model{T,D,Z,Y}
    prob::NLODEProblem{T,D,Z,Y}
    f::Function   
end =#
function save_prob_to_model(prob::NLODEProblem{T,D,Z,Y},path::String) where {T,D,Z,Y}
   # isdefined(Main,Symbol(modelName)) && error("model exists!")
#=     f=@RuntimeGeneratedFunction(prob.eqs)
    mymodel=Model{T,D,Z,Y}(prob,f) =#
    open(path, "a") do io
        # println(io,string(prob.eqs))
        
         println(io,string(prob.eqs))  
         

     end
end


function QSS_Solve_from_model(f::Function,prob::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V},initialTime::Float64,dQmin::Float64,dQrel::Float64,savetimeincrement::Float64) where {T,D,Z,Y,V}
    #f=prob.f
    #prob=model.prob
   # f=@RuntimeGeneratedFunction(prob.eqs)
    QSS_Unique_Solve(f,prob,finalTime,Val(V),initialTime,dQmin,dQrel,savetimeincrement)
end

function QSS_Solve(prob::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V},initialTime::Float64,dQmin::Float64,dQrel::Float64,savetimeincrement::Float64) where {T,D,Z,Y,V}
     f=@RuntimeGeneratedFunction(prob.eqs)
     QSS_Unique_Solve(f,prob,finalTime,Val(V),initialTime,dQmin,dQrel,savetimeincrement)
 end
 
 
function QSS_Unique_Solve(f::Function,prob::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V},initialTime::Float64,dQmin::Float64,dQrel::Float64,savetimeincrement::Float64) where {T,D,Z,Y,V}
    # ex=quote
     
     order=getOrderfromSolverMethod(Val(V))# this is to make a difference b/w qss1 and liqss1 but if performance then have qss_solve for liqss?
  
     # states = computeStates(prob.initConditions)
   # numberZC=size(prob.ZC_jacobian, 1)  #later test prob.Z
     quantum =  zeros(T)
     x = Vector{Taylor0{Float64}}(undef, T)
    # derx = Vector{Taylor0{Float64}}(undef, T)
     q = Vector{Taylor0{Float64}}(undef, T)
     nextStateTime = @MVector zeros(T)
     nextInputTime =  zeros(T)
     tx = @MVector zeros(T)
     tq = @MVector zeros(T)
     nextEventTime=@MVector zeros(Z)
     t = Taylor0(zeros(order + 1), order)
     t[1]=1.0
     integratorCache=Taylor0(zeros(order+1),order) #for integratestate only
     sizehint=floor(Int64, 1.0+(finalTime/savetimeincrement))
     savedVars = Vector{Array{Taylor0{Float64}}}(undef, T)# has to be vector (not SA) cuz to be resized in integrator
     #temparr=Array{Taylor0{Float64}}(undef, sizehint)
     for i = 1:T
     #push!(savedVars, zeros(20))# stiffness hint by                   user;& num_states
         temparr=Array{Taylor0{Float64}}(undef, sizehint)
         for j=1:sizehint          
             temparr[j]=Taylor0(zeros(order+1),order) # this number can be found from ft and saveat (ft/saveat=5/0.1=50) and maybe *2/3 factor (time can jump past a saving time) of user or default saveat if not provided
         end
         savedVars[i]=temparr
         x[i]=Taylor0(zeros(order + 1), order) 
         x[i][0]= prob.initConditions[i]
         # push!(q, Taylor0(zeros(order), order - 1))
         q[i]=Taylor0(zeros(order+1), order)#q normally 1order lower than x but since we want f(q) to  be a taylor that holds all info (1,2,3), lets have q of
        # derx[i]=Taylor0(zeros(order+1), order)
         #push!(savedVars[i], Taylor0(zeros(order + 1), order))
         
         tx[i] = initialTime
         tq[i] = initialTime
     end
     #savedTimes = Array{Float64}([initTime])
     #@show savedVars
     savedTimes = zeros(sizehint) #estimate like sizehint...later stiffness hint...to be multiplied by a stiffnessfactor
     savedTimes[1]=initialTime
 
 
     cacheSize=prob.cacheSize
     taylorOpsCache=Array{Taylor0{Float64},1}()
     for i=1:cacheSize
       push!(taylorOpsCache,Taylor0(zeros(order+1),order))
     end
 
    
 #=    zcf= @RuntimeGeneratedFunction(prob.zceqs) #args[2] cuz there is extra stuff
     eventf=@RuntimeGeneratedFunction(prob.eventEqus)  =#
 
     
     qssdata= QSS_data(quantum,x,q,tx,tq,nextStateTime,nextInputTime ,nextEventTime , t, integratorCache,order,savedVars,savedTimes,taylorOpsCache,finalTime,savetimeincrement, initialTime,dQmin,dQrel)
    
    if V==1 || V ==2 || V ==3
     QSS_integrate(Val(V),qssdata,prob,f)   # solver=Val(\d) needed to specialize the integrator
    else
     LiQSS_integrate(Val(V-3),qssdata,prob,f)
    end
     #return nothing be careful to add this
    # end
    # ex
 end
 
 
 function QSS_Solve(prob::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V}) where {T,D,Z,Y,V}
     initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1
     QSS_Solve(prob,finalTime,Val(V),initialTime,dQmin,dQrel,savetimeincrement)
 end
 function QSS_Solve(prob::NLODEProblem{T,D,Z,Y},finalTime::Float64) where {T,D,Z,Y}
     initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1
     QSS_Solve(prob,finalTime,Val(1),initialTime,dQmin,dQrel,savetimeincrement)
 end
 function QSS_Solve(prob::NLODEProblem{T,D,Z,Y}) where {T,D,Z,Y}
     initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1;finalTime=5.0
     QSS_Solve(prob,finalTime,Val(1),initialTime,dQmin,dQrel,savetimeincrement)
 end
 

 function QSS_Solve_from_model(f::Function,m::NLODEProblem{T,D,Z,Y},finalTime::Float64,::Val{V}) where {T,D,Z,Y,V}
    initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1
    QSS_Solve_from_model(f,m,finalTime,Val(V),initialTime,dQmin,dQrel,savetimeincrement)
end
function QSS_Solve_from_model(f::Function,m::NLODEProblem{T,D,Z,Y},finalTime::Float64) where {T,D,Z,Y}
    initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1
    QSS_Solve_from_model(f,m,finalTime,Val(1),initialTime,dQmin,dQrel,savetimeincrement)
end
function QSS_Solve_from_model(f::Function,m::NLODEProblem{T,D,Z,Y}) where {T,D,Z,Y}
    initialTime=0.0;dQmin=1e-6;dQrel=1e-3;savetimeincrement=0.1;finalTime=5.0
    QSS_Solve_from_model(f,m,finalTime,Val(1),initialTime,dQmin,dQrel,savetimeincrement)
end



