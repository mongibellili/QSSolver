
#helper struct that holds dependency metadata of an event (which vars exist on which side lhs=rhs

struct EventDependencyStruct
    id::Int
    evCont::Vector{Int}
    evDisc::Vector{Int}
    evContRHS::Vector{Int}
end

#struct that holds prob data 
struct NLODEDiscProblem{PRTYPE,T,Z,Y}<: NLODEProblem{PRTYPE,T,Z,Y} 
    prtype::Val{PRTYPE}
    a::Val{T}
    b::Val{Z}
    c::Val{Y}
    cacheSize::Int
    initConditions::Vector{Float64}  
    discreteVars::Vector{Float64}   
    jac::Vector{Vector{Int}}#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD
    ZCjac::Vector{Vector{Int}}
    eqs::Function#function that holds all ODEs
    eventDependencies::Vector{EventDependencyStruct}# 
    SD::Vector{Vector{Int}}#  I have a var and I want the der that are affected by it
    HZ::Vector{Vector{Int}}#  an ev occured and I want the ZC that are affected by it
    HD::Vector{Vector{Int}}#  an ev occured and I want the der that are affected by it
    SZ::Vector{Vector{Int}}#  I have a var and I want the ZC that are affected by it
    map::Function
end
struct savedNLODEDiscProblem{PRTYPE,T,Z,Y}<:  NLODEProblem{PRTYPE,T,Z,Y} 
    prtype::Val{PRTYPE}
    cacheSize::Int
    #initConditions::SVector{T,Float64}    
    initConditions::Vector{Float64}  
    discreteVars::Vector{Float64}   
    jacInts::Vector{Vector{Int}}#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD
   # eqs::Expr#function that holds all ODEs
    zc_SimpleJac::SVector{Z,SVector{T,Int}}#ZC jac dependency...later do not use SA for large sys...opposite of SZ
    #ZC_jacDiscrete::SVector{Z,SVector{D,Basic}}
    eventDependencies::SVector{Y,EventDependencyStruct}# 
    #initJac::MVector{T,MVector{T,Float64}}
    initJac :: Vector{Vector{Float64}} #initial values of jac
    #SD::SVector{T,SVector{T,Int}}
    SD::Vector{Vector{Int}}#  I have a var and I want the der that are affected by it
    HZ::SVector{Y,SVector{Z,Int}}#  an ev occured and I want the ZC that are affected by it
    HD::SVector{Y,SVector{T,Int}}#  an ev occured and I want the der that are affected by it
    SZ::SVector{T,SVector{Z,Int}}#  I have a var and I want the ZC that are affected by it
end

function getInitCond(prob::NLODEDiscProblem,i::Int)
    return prob.initConditions[i]
end
function getInitCond(prob::savedNLODEDiscProblem,i::Int)
    return prob.initConditions(i)
end

# receives user code and creates the problem struct
function NLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{Z},initCond::Vector{Float64},du::Symbol)where {T,D,Z}
    if verbose println("discrete nlodeprobfun  T D Z= $T $D $Z") end
    discrVars=Vector{Float64}()
    equs=Dict{Union{Int,Expr},Expr}()
    jac = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()# set used because do not want to insert an existing varNum
    SD = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()  # NO need to constrcut SD...SDVect will be extracted from Jac (needed for func-save case)
    jacDiscrete = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()# for now is similar to continous jac...if feature of d[i] not to be added then  jacDiscr=Dict{Union{Int,Expr},Set{Int}}()...later
    dD = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()  # similarly if feature not given then dD=Dict{Int,Set{Union{Int,Symbol,Expr}}}()

    zcequs=Vector{Expr}()#vect to collect if-statements
    eventequs=Vector{Expr}()#vect to collect events

    ZCjac=Vector{Vector{Int}}()#zeros(SVector{Z,SVector{T,Int}}) # bad when T is large...it is a good idea to use Vector{Vector{Int}}(undef, Z)
    ZCjacDiscrete=Vector{Vector{Int}}()#zeros(SVector{Z,SVector{D,Int}})# usually D (num of disc vars) and Z (num of zcfunctions) are small even for real problems
    ZCFCounter=0
    #SZ=Vector{Vector{Int}}()#zeros(SVector{T,SVector{Z,Int}}) # bad when T is large...(opposite of ZCjac)...many vars dont influence zcf so it is better to use dict
   # dZ=Vector{Vector{Int}}()#zeros(SVector{D,SVector{Z,Int}})# usually D (num of disc vars) and Z (num of zcfunctions) are small even for real problems (opposite of ZCjacDiscrete)
   SZ= Dict{Int,Set{Int}}()
   dZ= Dict{Int,Set{Int}}()  #since D is not large ...i can use zeros(SVector{D,SVector{Z,Int}})
   evsArr = EventDependencyStruct[]
    num_cache_equs=1#cachesize
    for argI in odeExprs.args
        #only diff eqs: du[]= number/one ref/call  
        if argI isa Expr &&  argI.head == :(=) && argI.args[1]== :discrete
            discrVars = Vector{Float64}(argI.args[2].args)
        elseif argI isa Expr &&  argI.head == :(=)  && argI.args[1] isa Expr && argI.args[1].head == :ref && argI.args[1].args[1]==du#&& ((argI.args[2] isa Expr && (argI.args[2].head ==:ref || argI.args[2].head ==:call ))||argI.args[2] isa Number)
            y=argI.args[1];rhs=argI.args[2]
            varNum=y.args[2] # order of variable
            if rhs isa Number # rhs of equ =number  
                equs[varNum]=:($((transformFSimplecase(:($(rhs))))))
            elseif rhs.head==:ref #rhs is only one var
               # if rhs.args[1]==:q # check not needed
                extractJacDepNormal(varNum,rhs,jac,jacDiscrete ,SD,dD ) 
               # end
                equs[varNum ]=:($((transformFSimplecase(:($(rhs))))))
            else #rhs head==call...to be tested later for  math functions and other possible scenarios or user erros                 
                extractJacDepNormal(varNum,rhs,jac,jacDiscrete ,SD,dD ) 
                temp=(transformF(:($(rhs),1))).args[2]  #number of caches distibuted   ...no need interpolation and wrap in expr....before was cuz quote....
                if num_cache_equs<temp 
                        num_cache_equs=temp
                end 
                equs[varNum]=rhs
            end 
        elseif @capture(argI, for counter_ in b_:niter_ loopbody__ end)
             specRHS=loopbody[1].args[2]
             extractJacDepLoop(b,niter,specRHS,jac ,jacDiscrete ,SD,dD  )   
             temp=(transformF(:($(specRHS),1))).args[2]
                if num_cache_equs<temp 
                    num_cache_equs=temp
                end 
               equs[:(($b,$niter))]=specRHS            
        elseif argI isa Expr && argI.head==:if   #@capture if did not work
           zcf=argI.args[1].args[2]
           ZCFCounter+=1
           extractZCJacDepNormal(ZCFCounter,zcf,ZCjac ,ZCjacDiscrete ,SZ ,dZ ) 
            if zcf.head==:ref  #if one_Var
                  push!(zcequs,(transformFSimplecase(:($(zcf)))))    
                ########################push!(zcequs,:($((transformFSimplecase(:($(zcf))))))) 
              #  push!(num_cache_zcequs,1) #to be deleted later 
            else # if whole expre ops with many vars            
                temp=:($((transformF(:($(zcf),1))).args[2]))   #number of caches distibuted, given 1 place holder for ex.args[2] to be filled inside and returned
                if num_cache_equs<temp 
                      num_cache_equs=temp
                end 
                push!(zcequs,(zcf))      
            end     
            ##################################################################################################################
            #                                                      events     
            ##################################################################################################################                  
            # each 'if-statmets' has 2 events (arg[2]=posEv and arg[3]=NegEv) each pos or neg event has a function...later i can try one event for zc
            if length(argI.args)==2  #if user only wrote the positive evnt, here I added the negative event wich does nothing
                nothingexpr = quote nothing end # neg dummy event 
                push!(argI.args, nothingexpr)
                Base.remove_linenums!(argI.args[3])
            end
      #  end#end temporary if
      # end#end temporay for argI
     
            #pos event
            newPosEventExprToFunc=changeExprToFirstValue(argI.args[2])  #change u[1] to u[1][0]  # pos ev can't be a symbol ...later maybe add check anyway
            push!(eventequs,newPosEventExprToFunc) 
            #neg eve
            if argI.args[3].args[1] isa Expr # argI.args[3] isa expr and is either an equation or :block symbol end ...so lets check .args[1]
                newNegEventExprToFunc=changeExprToFirstValue(argI.args[3])
                push!(eventequs,newNegEventExprToFunc) 
            else
                push!(eventequs,argI.args[3]) #symbol nothing
            end
            #after constructing the equations we move to dependencies: we need to change A[n] to An so that they become symbols
            posEvExp =  argI.args[2]
            negEvExp =  argI.args[3]
           #dump(negEvExp)
            indexPosEv = 2 * length(zcequs) - 1 # store events in order
            indexNegEv = 2 * length(zcequs)   
              #------------------pos Event--------------------#
            posEv_disArrLHS= Vector{Int}()#@SVector fill(NaN, D)   #...better than @SVector zeros(D), I can use NaN
            posEv_conArrLHS= Vector{Int}()#@SVector fill(NaN, T)  
            posEv_conArrRHS=Vector{Int}()#@SVector zeros(T)    #to be used inside intgrator to updateOtherQs (intgrateState) before executing the event there is no discArrRHS because d is not changing overtime to be updated      
            for j = 1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
                !(posEvExp.args[j]  isa Expr &&  posEvExp.args[j].head == :(=))  && error("event should be A=B")
                poslhs=posEvExp.args[j].args[1];posrhs=posEvExp.args[j].args[2]
                !(poslhs  isa Expr &&  poslhs.head == :ref && (poslhs.args[1]==:q || poslhs.args[1]==:d)) && error("lhs of events must be a continuous or a discrete variable")
                if poslhs.args[1]==:q
                    push!(posEv_conArrLHS,poslhs.args[2])
                else # lhs is a disc var 
                    push!(posEv_disArrLHS,poslhs.args[2])
                end
                #@show poslhs,posrhs
                postwalk(posrhs) do a   #
                    if a isa Expr && a.head == :ref && a.args[1]==:q# 
                        push!(posEv_conArrRHS,  (a.args[2]))  #                    
                    end
                    return a 
                end

            end
            #------------------neg Event--------------------#
            negEv_disArrLHS= Vector{Int}()#@SVector fill(NaN, D)   #...better than @SVector zeros(D), I can use NaN
            negEv_conArrLHS= Vector{Int}()#@SVector fill(NaN, T)  
            negEv_conArrRHS=Vector{Int}()#@SVector zeros(T)    #to be used inside intgrator to updateOtherQs (intgrateState) before executing the event there is no discArrRHS because d is not changing overtime to be updated      
            if negEvExp.args[1] != :nothing
                for j = 1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
                    !(negEvExp.args[j].args  isa Expr &&  argI.head == :(=))  && error("event should be A=B")
                    neglhs=negEvExp.args[j].args[1];negrhs=negEvExp.args[j].args[1]
                    !(neglhs  isa Expr &&  neglhs.head == :ref && (neglhs.args[1]==:q || neglhs.args[1]==:q)) && error("lhs of events must be a continuous or a discrete variable")
                    if neglhs.args[1]==:q
                        push!(negEv_conArrLHS,neglhs.args[2])
                    else # lhs is a disc var 
                        push!(negEv_disArrLHS,neglhs.args[2])
                    end
                    negtwalk(negrhs) do a   #
                        if a isa Expr && a.head == :ref && a.args[1]==:q# 
                            push!(negEv_conArrRHS,  (a.args[2]))  #                    
                        end
                        return a 
                    end
            
                end
            end 
            structposEvent = EventDependencyStruct(indexPosEv, posEv_conArrLHS, posEv_disArrLHS,posEv_conArrRHS) # right now posEv_conArr is vect of floats
            push!(evsArr, structposEvent)
            structnegEvent = EventDependencyStruct(indexNegEv, negEv_conArrLHS, negEv_disArrLHS,negEv_conArrRHS)
            push!(evsArr, structnegEvent)
        end #end cases 
 
    end #end for #########################################################################################################
   ############@show evsArr
  
    #println("end of useless convesrion to svector")
    allEpxpr=Expr(:block)
    ##############diffEqua######################
   #=  io = IOBuffer() # i guess this is a sys of diff equ solver so I am not optimizing for case 1 equ
    write(io, "if j==1  $(equs[1]) ;return nothing") =#
    s="if i==0 return nothing\n"  # :i is the mute var
    for elmt in equs
        Base.remove_linenums!(elmt[1])
        Base.remove_linenums!(elmt[2])
        if elmt[1] isa Int
            s*="elseif i==$(elmt[1]) $(elmt[2]) ;return nothing\n"
        end
        if elmt[1] isa Expr
            s*="elseif $(elmt[1].args[1])<=i<=$(elmt[1].args[2]) $(elmt[2]) ;return nothing\n"
        end
    end
    s*=" end "
  
    myex1=Meta.parse(s)
    push!(allEpxpr.args,myex1)

     ##############ZCF######################
     
    if length(zcequs)>0
       #=  io = IOBuffer()
        write(io, "if zc==1  $(zcequs[1]) ;return nothing") =#
        s="if zc==1  $(zcequs[1]) ;return nothing"
        for i=2:length(zcequs)
            s*= " elseif zc==$i $(zcequs[i]) ;return nothing"
        end
        # write(io, " else  $(zcequs[length(zcequs)]) ;return nothing end ")
        s*= " end "
       
        myex2=Meta.parse(s)
        push!(allEpxpr.args,myex2)
    end
    ##############events######################
    if length(eventequs)>0
       
        s= "if ev==1  $(eventequs[1]) ;return nothing"
        for i=2:length(eventequs)
            s*= " elseif ev==$i $(eventequs[i]) ;return nothing"
        end
        # write(io, " else  $(zcequs[length(zcequs)]) ;return nothing end ")
        s*= " end "
       
        myex3=Meta.parse(s)
        push!(allEpxpr.args,myex3)
    end

    Base.remove_linenums!(allEpxpr)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :f   
    def[:args] = [:(i::Int),:(zc::Int),:(ev::Int),:(q::Vector{Taylor0}),:(d::Vector{Float64}), :(t::Taylor0),:(cache::Vector{Taylor0})]
    def[:body] = allEpxpr
    #def[:rtype]=:nothing# test if cache1 always holds sol  
    functioncode=combinedef(def)
    
  


#=   
jac = Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set([2]))     
SD = Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(2 => Set([1]))      
jacDiscrete = Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}(1 => Set())
dD = Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}}()
ZCjac = [[1], [2]]
ZCjacDiscrete = [[1], Int64[]]
ZCFCounter = 2
SZ = Dict{Int64, Set{Int64}}(2 => Set([2]), 1 => Set([1]))
dZ = Dict{Int64, Set{Int64}}(1 => Set([1])) =#

jacVect=createJacVect(jac,Val(T))
SDVect=createSDVect(jac,Val(T))
dDVect =createdDvect(dD)
jacDiscreteVect=createJacVect(jacDiscrete,Val(T))
SZvect=createSZvect(SZ,Val(T))

#= SDVect = [Int64[], [2, 1], [3], [4], Int64[]]    
jacVect = [[2], [2], [3], [4], Int64[]]
jacDiscreteVect = [[2, 1], [1], [1], [1], [2, 1]]
dDVect = [[5, 2, 3, 4, 1], [5, 1]] =#

    # temporary dependencies to be used to determine HD and HZ...determine HD: event-->Derivative   && determine HZ:Event-->ZCfunction....influence of events on derivatives and zcfunctions:
    #an event is a discrteVar change or a cont Var change. So HD=HD1 UNION HD2  (same for HZ=HZ1 UNION HZ2)
    # (1) through a discrete Var: 
    #  ============================
    # HD1=Hd-->dD  where Hd  comes from the eventDependecies Struct. 

    # HZ1=Hd-->dZ  where Hd  comes from the eventDependecies Struct. 
    HZ1HD1=createDependencyToEventsDiscr(dDVect,dZ,evsArr) 

    # (2) through a continous Var: 
    # ==============================
     #  HD2=Hs-->sD where  Hs comes from the eventDependecies Struct.(sd already created) 
     #  HZ2=Hs-->sZ where  Hs comes from the eventDependecies Struct.(sZ already created) 
    HZ2HD2=createDependencyToEventsCont(SDVect,SZ,evsArr) 
   # @show HZ2HD2
    ##########UNION##############
    HZ=unionDependency(HZ1HD1[1],HZ2HD2[1])
    HD=unionDependency(HZ1HD1[2],HZ2HD2[2])

 
    mapFun=createMapFun(jac,:f)
    mapFunF=@RuntimeGeneratedFunction(mapFun)
    functioncodeF=@RuntimeGeneratedFunction(functioncode)
     
    
    myodeProblem = NLODEDiscProblem(Val(1),Val(T),Val(Z),Val(2Z),num_cache_equs,initCond, discrVars, jacVect ,ZCjac  ,functioncodeF, evsArr,SDVect,HZ,HD,SZvect,mapFunF)
  

end


function saveNLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{D},::Val{Z},initCond::Dict{Union{Int,Expr},Int},du::Symbol)where {T,D,Z}
    prob=NLodeProblemFunc(odeExprs,Val(T),Val(D),Val(Z))
    def=splitdef(prob.eqs)
    if (odeExprs.args[1] isa Symbol) 
        #isdefined(Main,odeExprs.args[1]) && error("model exists!")
        def[:name]=odeExprs.args[1]
    end   
    functioncode=combinedef(def)
    allEpxpr=Expr(:block)

    s= " zc_SimpleJac=SVector{$Z,SVector{$T,Int}}(tuple($(prob.zc_SimpleJac)...))"
    myex1=Meta.parse(s)
    push!(allEpxpr.args,myex1)

    s= " eventDependencies=SVector{$(2Z),EventDependencyStruct}(tuple($(prob.eventDependencies)...))"
    myex1=Meta.parse(s)
    push!(allEpxpr.args,myex1)

    s= " HZ=SVector{$(2Z),SVector{$Z,Int}}(tuple($(prob.HZ)...))"
    myex1=Meta.parse(s)
    push!(allEpxpr.args,myex1)

    s= " HD=SVector{$(2Z),SVector{$T,Int}}(tuple($(prob.HD)...))"
    myex1=Meta.parse(s)
    push!(allEpxpr.args,myex1)

    s= " SZ=SVector{$(T),SVector{$Z,Int}}(tuple($(prob.SZ)...))"
    myex1=Meta.parse(s)
    push!(allEpxpr.args,myex1)

    s="myodeProblem = savedNLODEDiscProblem($(prob.cacheSize),$(prob.initConditions), $(prob.discreteVars), $(prob.jacInts), zc_SimpleJac, eventDependencies,$(prob.initJac),$(prob.SD),HZ,HD,SZ)"

   


    myex1=Meta.parse(s)
    push!(allEpxpr.args,myex1)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    if !(odeExprs.args[1] isa Symbol) 
        def[:name] = :fproblem   
    else
        #isdefined(Main,odeExprs.args[1]) && error("model exists!")
        def[:name]=Symbol(odeExprs.args[1],:_prob)
    end
    def[:body] = allEpxpr
    problemcode=combinedef(def)
    open(odeExprs.args[2], "a") do io        
        println(io,string(functioncode))  
        println(io,string(problemcode)) 
    end
    return prob # just in case you want to solve/run while saving a problem    
end


function createdDvect(dD::Dict{Union{Int64, Expr}, Set{Union{Int64, Expr, Symbol}}})
    lendD = length(dD)
    dDVect = Vector{Vector{Int}}(undef, lendD)
    for ii=1:lendD
        dDVect[ii]=Vector{Int}()# define it so i can push elements as i find them below
    end
    for dictElement in dD
        temp=Vector{Int}()
        for k in dictElement[2] 
            if k isa Int
                push!(temp,k)
            elseif k isa Expr #tuple
                for j_=(k.args[1]):(k.args[2])  
                    push!(temp,j_)
                end
            end
        end
       dDVect[dictElement[1]]=temp
    end
    dDVect
end


function createSZvect(SZ :: Dict{Int64, Set{Int64}},::Val{T})where{T}
    szVect = Vector{Vector{Int}}(undef, T)
    for ii=1:T
        szVect[ii]=Vector{Int}()# define it so i can push elements as i find them below
    end
    for ele in SZ
        szVect[ele[1]]=collect(ele[2])
    end
    szVect
end