

# user does not provide solver. default mliqss2
function QSS_Solve(prob::NLODEProblem{PRTYPE,T,Z,Y,CS};sparsity::Val{Sparsity}=Val(false),finalTime=1.0::Float64,saveat=0.1::Float64,initialTime=0.0::Float64,dQmin=1e-6::Float64,dQrel=1e-3::Float64,maxErr=Inf::Float64) where {PRTYPE,T,Z,Y,CS,Sparsity}    
   QSS_Solve(prob,(Val(:mliqss),Val(2));sparsity=sparsity,finalTime=finalTime,saveat=saveat,initialTime=initialTime,dQmin=dQmin,dQrel=dQrel,maxErr=maxErr)  
end
#main solve interface
function QSS_Solve(prob::NLODEProblem{PRTYPE,T,Z,Y,CS},::Tuple{SolverType, OrderType};sparsity::Val{Sparsity}=Val(false),finalTime=1.0::Float64,saveat=0.1::Float64,initialTime=0.0::Float64,dQmin=1e-6::Float64,dQrel=1e-3::Float64,maxErr=Inf::Float64) where{PRTYPE,T,Z,Y,CS,SolverType,OrderType,Sparsity}    
   custom_Solve(prob,SolverType,OrderType,Val(Sparsity),finalTime,saveat,initialTime,dQmin,dQrel,maxErr)
end
#default solve method: this is not to be touched...extension or modification is done through creating another custom_solve
function custom_Solve(prob::NLODEProblem{PRTYPE,T,Z,Y,CS},::Type{Val{Solver}},::Type{Val{Order}},::Val{Sparsity},finalTime::Float64,saveat::Float64,initialTime::Float64,dQmin::Float64,dQrel::Float64,maxErr::Float64) where{PRTYPE,T,Z,Y,CS,Solver,Order,Sparsity}
    sizehint=floor(Int64, 1.0+(finalTime/saveat)*0.6)
    commonQSSdata=createCommonData(prob,Val(T),Val(Z),Val(Order),sizehint,finalTime,saveat, initialTime,dQmin,dQrel,maxErr)
    jac=getClosure(prob.jac)::Function
    SD=getClosure(prob.SD)::Function
    if Solver==:qss
        QSS_integrate(commonQSSdata,prob,prob.eqs,jac,SD,prob.map)
    else
          liqssdata=createLiqssData(prob,Val(Sparsity),Val(T),Val(Order))
         specialLiqssData=createSpecialLiqssData(Val(T))
        if Solver==:nmliqss
             nmLiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.map)
        elseif Solver==:nliqss
            nLiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.map)
        elseif Solver==:mliqss
            mLiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.map)
        elseif Solver==:liqss
             LiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.map)     
        end
    end
 end



 function getClosure(jacSD::Function)::Function # if not used with saving, it allocates...happens when multiple fun_'combinedef' and runtimegenerated
   function closureJacSD(i::Int)
        jacSD(i)
   end
   return closureJacSD
 end

 function getClosure(jacSD::Vector{Vector{Int}})::Function
    function closureJacSD(i::Int)
         jacSD[i]
    end
    return closureJacSD
  end




#helper methods...not to be touched...extension can be done through creating others via specializing on one PRTYPE or more of the symbols (PRTYPE,T,Z,Y,Order) if in the future...
#################################################################################################################################################################################
function createCommonData(prob::NLODEProblem{PRTYPE,T,Z,Y,CS},::Val{T},::Val{Z},::Val{Order},sizehint::Int,finalTime::Float64,saveat::Float64,initialTime::Float64,dQmin::Float64,dQrel::Float64,maxErr::Float64)where{PRTYPE,T,Z,Y,CS,Order}
    quantum =  zeros(T)
    x = Vector{Taylor0}(undef, T)
    q = Vector{Taylor0}(undef, T)
    savedTimes=Vector{Vector{Float64}}(undef, T)
    savedVars = Vector{Vector{Float64}}(undef, T)
    nextStateTime =  zeros(T)
    nextInputTime =   zeros(T)
    tx = zeros(T)
    tq =  zeros(T)
    nextEventTime=@MVector zeros(Z)# only Z number of zcf is usually under 100...so use of SA is ok
    
    t = Taylor0(zeros(Order + 1), Order)
    t[1]=1.0
    t[0]=initialTime
    integratorCache=Taylor0(zeros(Order+1),Order) #for integratestate only

    for i = 1:T
        nextInputTime[i]=Inf
        x[i]=Taylor0(zeros(Order + 1), Order) 
        x[i][0]= getInitCond(prob,i)        # x[i][0]= prob.initConditions[i] if to remove saving as func
        q[i]=Taylor0(zeros(Order+1), Order)#q normally 1order lower than x but since we want f(q) to  be a taylor that holds all info (1,2,3), lets have q of same Order and not update last coeff        
        tx[i] = initialTime
        tq[i] = initialTime
        savedTimes[i]=Vector{Float64}()
        savedVars[i]=Vector{Float64}()
    end
    
    taylorOpsCache=Array{Taylor0,1}()# cache= vector of taylor0s of size CS
    for i=1:CS
    push!(taylorOpsCache,Taylor0(zeros(Order+1),Order))
    end
    
    commonQSSdata= CommonQSS_data(Val(Order),quantum,x,q,tx,tq,nextStateTime,nextInputTime ,nextEventTime , t, integratorCache,taylorOpsCache,finalTime,saveat, initialTime,dQmin,dQrel,maxErr,savedTimes,savedVars)
end





#no sparsity
function createLiqssData(prob::NLODEProblem{PRTYPE,T,Z,Y,CS},::Val{false},::Val{T},::Val{Order})where{PRTYPE,T,Z,Y,CS,Order}
    a = Vector{Vector{Float64}}(undef, T)
   # u=Vector{Vector{MVector{Order,Float64}}}(undef, T)
    qaux = Vector{MVector{Order,Float64}}(undef, T)
    dxaux=Vector{MVector{Order,Float64}}(undef, T)
    olddx = Vector{MVector{Order,Float64}}(undef, T)
    olddxSpec = Vector{MVector{Order,Float64}}(undef, T)
    #= @timeit  "liqssdense" =# for i=1:T
       #=  temparr=Vector{MVector{Order,Float64}}(undef, T)
        for j=1:T
            temparr[j]=zeros(MVector{Order,Float64})
        end
        u[i]=temparr =#
        a[i]=zeros(T)
        qaux[i]=zeros(MVector{Order,Float64})
        olddx[i]=zeros(MVector{Order,Float64})
        dxaux[i]=zeros(MVector{Order,Float64})
        olddxSpec[i]=zeros(MVector{Order,Float64})

    end
    liqssdata= LiQSS_data(Val(false),a#= ,u =#,qaux,olddx,dxaux,olddxSpec)
end

#to be removed if sparsity did not help
function createLiqssData(prob::NLODEProblem{PRTYPE,T,Z,Y,CS},::Val{true},::Val{T},::Val{Order})where{PRTYPE,T,Z,Y,CS,Order}
    a = Vector{Vector{Float64}}(undef, T)
    u=Vector{Vector{MVector{Order,Float64}}}(undef, T)
    qaux = Vector{MVector{Order,Float64}}(undef, T)
    olddx = Vector{MVector{Order,Float64}}(undef, T)
    olddxSpec = Vector{MVector{Order,Float64}}(undef, T)
    jacDim=prob.jacDim
    #= @timeit  "liqsssparse" =# for i=1:T
        temparr=Vector{MVector{Order,Float64}}(undef, T)
        for j=1:T
            temparr[j]=zeros(MVector{Order,Float64})
        end
        u[i]=temparr
        a[i]=zeros(jacDim(i))
        qaux[i]=zeros(MVector{Order,Float64})
        olddx[i]=zeros(MVector{Order,Float64})
        olddxSpec[i]=zeros(MVector{Order,Float64})

    end
    liqssdata= LiQSS_data(Val(true),a,u,qaux,olddx,olddxSpec)
end


function createSpecialLiqssData(::Val{T})where{T}
    cacheA=zeros(MVector{1,Int})
    direction= zeros(T)
    qminus= zeros(T)
    buddySimul=zeros(MVector{2,Int})
    prevStepVal = zeros(T)
    specialLiqssData=SpecialLiQSS_data(cacheA,direction,qminus,buddySimul,prevStepVal)
end





# get init conds for normal vect...getinitcond for fun can be found with qssnlsavedprob file
function getInitCond(prob::NLODEContProblem,i::Int)
    return prob.initConditions[i]
end





















