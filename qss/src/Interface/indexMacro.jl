

macro NLodeProblem(odeExprs)
    Base.remove_linenums!(odeExprs)
    if verbose println("starting prob parsing...") end 
    probHelper=arrangeProb(odeExprs)# replace symbols and params , extract info about size,symbols,initconds
    probSize=probHelper.problemSize
    discSize=probHelper.discreteSize
    zcSize=probHelper.numZC
    initConds=probHelper.initConditions # vector
    du=probHelper.du
    if length(initConds)==0  #user chose shortcuts...initcond saved in a dict
        initConds=zeros(probSize)# vector of init conds to be created
        savedinitConds=probHelper.savedInitCond #init conds already created in a dict
        for i in savedinitConds
            if i[1] isa Int #if key (i[1]) is an integer
            initConds[i[1]]=i[2]
            else # key is an expression
                for j=(i[1].args[1]):(i[1].args[2])  
                    initConds[j]=i[2]
                end
            end
        end
    end

    NLodeProblemFunc(odeExprs,Val(probSize),Val(discSize),Val(zcSize),initConds,du)     #returns continuours prob   
end

#= macro saveNLodeProblem(odeExprs) # recommended for large cont problems to save first
    Base.remove_linenums!(odeExprs)
    if verbose println("starting prob saving...") end 
    probHelper=arrangeProb(odeExprs)
    probSize=probHelper.problemSize
    discSize=probHelper.discreteSize
    zcSize=probHelper.numZC
    initConds=probHelper.savedInitCond
    du=probHelper.du
    if length(initConds)==0 # user chose to represent initCond in a large vector eventhough savedNLode is for large prob and representing initconds in a dict is more efficient
        vectInitCond=probHelper.initConditions
        for i=1:length(vectInitCond)
            initConds[i]=vectInitCond[i]
        end
    end
    saveNLodeProblemFunc(odeExprs,Val(probSize),Val(discSize),Val(zcSize),initConds,du) 
end
 =#
struct probHelper #helper struct to return stuff from arrangeprob
    problemSize::Int
    discreteSize::Int
    numZC::Int
    savedInitCond::Dict{Union{Int,Expr},Float64}
    initConditions::Vector{Float64}
    du::Symbol
end
function arrangeProb(x::Expr) # replace symbols and params , extract info about size,symbols,initconds
    param=Dict{Symbol,Number}()
    stateVarName=:q
    du=:nothing #default anything 
    problemSize=0
    discreteSize=0
    numZC=0
    savedInitCond=Dict{Union{Int,Expr},Float64}()
    initConditions=Vector{Float64}()
    for argI in x.args
        if argI isa Expr &&  argI.head == :(=) #&& @capture(argI, y_ = rhs_) 
            y=argI.args[1];rhs=argI.args[2]
            if y isa Symbol && rhs isa Number  #params: fill dict of param
                param[y]=rhs
            elseif y isa Expr && y.head == :ref && rhs isa Number #initial conds "1st way"  u[a]=N or u[a:b]=N...
               if string(y.args[1])[1] !='d' #prevent the case diffEq du[]=Number
                    stateVarName=y.args[1]   #extract var name
                    if du==:nothing du=Symbol(:d,stateVarName) end # 
                    if  y.args[2] isa Expr && y.args[2].args[1]== :(:)  #u[a:b]=N
                            y.args[2].args!=3 || error(" use syntax u[a:b]") # not needed
                            savedInitCond[:(($(y.args[2].args[2]),$(y.args[2].args[3])))]=rhs# dict {expr->float}
                            if problemSize < y.args[2].args[3]
                                problemSize=y.args[2].args[3]
                            end
                    elseif y.args[2] isa Int   #u[a]=N
                        problemSize+=1
                        savedInitCond[y.args[2]]=rhs#.args[1]
                    end 
                end  
            elseif y isa Symbol && rhs isa Expr && rhs.head==:vect # cont vars u=[] "2nd way" or discrete vars d=[]
                if y!=:discrete
                    stateVarName=y
                    du=Symbol(:d,stateVarName)
                    initConditions= convert(Array{Float64,1}, rhs.args)
                    problemSize=length(rhs.args)
                else
                    discreteSize = length(rhs.args)  
                end    
            elseif y isa Expr && y.head == :ref && rhs isa Expr && rhs.head !=:vect#&& rhs.head==:call # a diff equa not in a loop
                argI.args[2]=changeVarNames_params(rhs,stateVarName,:nothing,param)
            end
        elseif @capture(argI, for var_ in b_:niter_ loopbody__ end)
            argI.args[2]=changeVarNames_params(loopbody[1],stateVarName,var,param)
        elseif argI isa Expr && argI.head==:if
            numZC+=1
            (length(argI.args)!=3 && length(argI.args)!=2) && error("use format if A>0 B else C or if A>0 B")
            !(argI.args[1] isa Expr && argI.args[1].head==:call && argI.args[1].args[1]==:> && (argI.args[1].args[3]==0||argI.args[1].args[3]==0.0)) && error("use the format 'if a>0: change if a>b to if a-b>0")
              !(argI.args[1].args[2] isa Expr) && error("LHS of >  must be be an expression!")
              argI.args[1].args[2]=changeVarNames_params(argI.args[1].args[2],stateVarName,:nothing,param)# name changes have to be per block
              if length(argI.args)==2 #user used if a b
                argI.args[2]=changeVarNames_params(argI.args[2],stateVarName,:nothing,param)
              elseif length(argI.args)==3 #user used if a b else c
                argI.args[2]=changeVarNames_params(argI.args[2],stateVarName,:nothing,param)
                argI.args[3]=changeVarNames_params(argI.args[3],stateVarName,:nothing,param)
                
              end
        end#end cases of argI
    end#end for argI in args
    p=probHelper(problemSize,discreteSize,numZC,savedInitCond,initConditions,du)
end#end function




qss1()=(Val(:qss),Val(1))
qss2()=(Val(:qss),Val(2))
qss3()=(Val(:qss),Val(3))

nmliqss1()=(Val(:nmliqss),Val(1))
nmliqss2()=(Val(:nmliqss),Val(2))
nmliqss3()=(Val(:nmliqss),Val(3))

nliqss1()=(Val(:nliqss),Val(1))
nliqss2()=(Val(:nliqss),Val(2))
nliqss3()=(Val(:nliqss),Val(3))

mliqss1()=(Val(:mliqss),Val(1))
mliqss2()=(Val(:mliqss),Val(2))
mliqss3()=(Val(:mliqss),Val(3))

liqss1()=(Val(:liqss),Val(1))
liqss2()=(Val(:liqss),Val(2))
liqss3()=(Val(:liqss),Val(3))



sparse()=Val(true)
dense()=Val(false)

