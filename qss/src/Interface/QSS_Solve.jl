

# user does not provide solver. default mliqss2
function QSS_Solve(prob::NLODEProblem{PRTYPE,T,Z,Y,CS};detection::Val{Detection}=Val(0),finalTime=1.0::Float64,saveat=0.1::Float64,initialTime=0.0::Float64,dQmin=1e-6::Float64,dQrel=1e-3::Float64,maxErr=Inf::Float64) where {PRTYPE,T,Z,Y,CS,Detection}    
    QSS_Solve(prob,(Val(:SimulAna),Val(2));detection=detection,finalTime=finalTime,saveat=saveat,initialTime=initialTime,dQmin=dQmin,dQrel=dQrel,maxErr=maxErr)  
 end
 #main solve interface
 function QSS_Solve(prob::NLODEProblem{PRTYPE,T,Z,Y,CS},::Tuple{SolverName, SolverType};detection::Val{Detection}=Val(0),finalTime=1.0::Float64,saveat=0.1::Float64,initialTime=0.0::Float64,dQmin=1e-6::Float64,dQrel=1e-3::Float64,maxErr=Inf::Float64) where{PRTYPE,T,Z,Y,CS,SolverName,SolverType,Detection}    
    custom_Solve(prob,SolverName,SolverType,Val(Detection),finalTime,saveat,initialTime,dQmin,dQrel,maxErr)
 end
#default solve method: this is not to be touched...extension or modification is done through creating another custom_solve with different PRTYPE
function custom_Solve(prob::NLODEProblem{PRTYPE,T,Z,Y,CS},::Type{Val{SolverName}},::Type{Val{SolverType}},::Val{Detection},finalTime::Float64,saveat::Float64,initialTime::Float64,dQmin::Float64,dQrel::Float64,maxErr::Float64) where{PRTYPE,T,Z,Y,CS,SolverName,SolverType,Detection}
    # sizehint=floor(Int64, 1.0+(finalTime/saveat)*0.6)
    commonQSSdata=createCommonData(prob,Val(T),Val(Z),Val(SolverName),Val(SolverType),finalTime,saveat, initialTime,dQmin,dQrel,maxErr)
    jac=getClosure(prob.jac)::Function #if in future jac and SD are different datastructures
    SD=getClosure(prob.SD)::Function
   
        liqssdata=createLiqssData(prob,Val(Detection),Val(T))
         specialLiqssData=createSpecialLiqssData(Val(T))
         integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)   
        #= if Solver==:nmliqss
             nmLiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)
        elseif Solver==:nliqss
            nLiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)
        elseif Solver==:mliqss
            mLiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)
        elseif Solver==:liqss
             LiQSS_integrate(commonQSSdata,liqssdata,specialLiqssData,prob,prob.eqs,jac,SD,prob.exactJac)     
        end =#
    
 end



 function getClosure(jacSD::Function)::Function # 
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




#helper methods...not to be touched...extension can be done through creating others via specializing on one PRTYPE or more of the symbols (PRTYPE,T,Z,D,Order) if in the future...
#################################################################################################################################################################################
function createCommonData(prob::NLODEProblem{PRTYPE,T,Z,Y,CS},::Val{T},::Val{Z},::Val{SolverName},::Val{SolverType},finalTime::Float64,saveat::Float64,initialTime::Float64,dQmin::Float64,dQrel::Float64,maxErr::Float64)where{PRTYPE,T,Z,Y,CS,SolverName,SolverType}
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
    
    t = Taylor0(zeros(1 + 1), 1)
    t[1]=1.0
    t[0]=initialTime
    integratorCache=Taylor0(zeros(1+1),1) #for integratestate only

    for i = 1:T
        nextInputTime[i]=Inf
        x[i]=Taylor0(zeros(1 + 1), 1) 
        x[i][0]= getInitCond(prob,i)        # x[i][0]= prob.initConditions[i] if to remove saving as func
        q[i]=Taylor0(zeros(1+1), 1)#q normally 1order lower than x but since we want f(q) to  be a taylor that holds all info (1,2,3), lets have q of same 1 and not update last coeff        
        tx[i] = initialTime
        tq[i] = initialTime
        savedTimes[i]=Vector{Float64}()
        savedVars[i]=Vector{Float64}()
    end
    
    taylorOpsCache=Array{Taylor0,1}()# cache= vector of taylor0s of size CS
    for i=1:CS
    push!(taylorOpsCache,Taylor0(zeros(1+1),1))
    end
    
    commonQSSdata= CommonQSS_data(Val(SolverName),Val(SolverType),quantum,x,q,tx,tq,nextStateTime,nextInputTime ,nextEventTime , t, integratorCache,taylorOpsCache,finalTime,saveat, initialTime,dQmin,dQrel,maxErr,savedTimes,savedVars)
end






function createLiqssData(prob::NLODEProblem{PRTYPE,T,Z,Y,CS},::Val{Detection},::Val{T})where{PRTYPE,T,Z,Y,CS,Detection}
    a = Vector{Vector{Float64}}(undef, T)
   # u=Vector{Vector{MVector{1,Float64}}}(undef, T)
    qaux = Vector{MVector{1,Float64}}(undef, T)
    dxaux=Vector{MVector{1,Float64}}(undef, T)
    olddx = Vector{MVector{1,Float64}}(undef, T)
    olddxSpec = Vector{MVector{1,Float64}}(undef, T)
    #= @timeit  "liqssdense" =# for i=1:T
       #=  temparr=Vector{MVector{1,Float64}}(undef, T)
        for j=1:T
            temparr[j]=zeros(MVector{1,Float64})
        end
        u[i]=temparr =#
        a[i]=zeros(T)
        qaux[i]=zeros(MVector{1,Float64})
        olddx[i]=zeros(MVector{1,Float64})
        dxaux[i]=zeros(MVector{1,Float64})
        olddxSpec[i]=zeros(MVector{1,Float64})

    end
    liqssdata= LiQSS_data(Val(Detection),a#= ,u =#,qaux,olddx,dxaux,olddxSpec)
end





function createSpecialLiqssData(::Val{T})where{T}
    cacheA=zeros(MVector{1,Float64})
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





















