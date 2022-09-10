
#helper struct that holds dependency data of an event
struct EventDependencyStruct{T,D} #<: AbstractEventDependecy
    id::Int
    evCont::SVector{T,Float64}
    evDisc::SVector{D,Float64}
end


struct NLODEProblem{T,D,Z,Y}
    initConditions::SVector{T,Float64}
    discreteVars::MVector{D,Float64}
    # jacobian::SVector{T,SVector{T,Float64}} 
    jacobian::SVector{T,SVector{T,Basic}}
    eqs::Vector{Expr}
   # dereqs::Vector{Expr}
   # derdereqs::Vector{Expr}
    zceqs::Vector{Expr}
    eventEqus::Vector{Expr}
    discreteJacobian::SVector{T,SVector{D,Basic}}
    #inputVars::SVector{T,Float64} 
    ZC_jacobian::SVector{Z,SVector{T,Basic}}
    ZC_jacDiscrete::SVector{Z,SVector{D,Basic}}
    # ZCinputVars::SVector{Z,Float64}   # in case input signal are function of t...use SVector{Z,Function} and create a different struct 
    eventDependencies::SVector{Y,EventDependencyStruct}# i think it can better if use matrix instead of EventDependencyStruct
end





macro NLodeProblem(odeExprs)
    Base.remove_linenums!(odeExprs)
    M = length(odeExprs.args)
    #check user input and throw error if not vector
    T = length((odeExprs.args[1].args[2]).args)  # number of cont vars.....correponds to the number of diff equ 
  # display((odeExprs.args[1].args[2]).args[1]/2)
    #check user input and throw error if not vector
    D = length((odeExprs.args[2].args[2]).args)   # number of discrete vars

    #check user input and throw error if not equation
    equs = [
        quote
            $(Symbol(:f_, i))(q::Vector{Taylor1{Float64}},d::Vector{Float64}, t::Taylor1{Float64}) = $(odeExprs.args[i].args[2])#anonymous func
        end for i = 3:2+T
    ]  #map(f ,l) is an alternative
    #zc = odeExpr.args[i].args[1].args[2]
    #check user input and throw error if not "if statement"
    ZCequs = [
        quote
            $(Symbol(:g_, i))(x::Vector{Taylor1{Float64}},d::Vector{Float64}, t::Taylor1{Float64}) = $(odeExprs.args[i].args[1].args[2])#anonymous func
        end for i = 3+T:M
    ]
   # println("equs= ",equs)
    eventEqus=Vector{Expr}()
    for i = 3+T:M# normally as a design decision these funcs should not have args...but i will keep them for now. also i probably should add return nothing
        push!(eventEqus,:($(Symbol(:h_, i))(q::Vector{Taylor1{Float64}},d::Vector{Float64}, t::Taylor1{Float64}) = $(odeExprs.args[i].args[2])))
        push!(eventEqus,:($(Symbol(:i_, i))(q::Vector{Taylor1{Float64}},d::Vector{Float64}, t::Taylor1{Float64}) = $(odeExprs.args[i].args[3])))
        
    end
#println("eventequs= ",eventEqus)
   # odeExpr = postwalk(x -> x isa Expr && x.head == :ref ? dump(x) : x, odeExprs)
    odeExpr = postwalk(x -> x isa Expr && x.head == :ref ? Symbol((x.args[1]), (x.args[2])) : x, odeExprs)
     #println("-----modified-----=",odeExpr)

    #--------------------------------------------------------
    #               intial vals and diff equa part           #
    #---------------------------------------------------------

    contVars = SVector{T,Float64}((odeExpr.args[1].args[2]).args)
    discrVars = MVector{D,Float64}((odeExpr.args[2].args[2]).args)

    jac = zeros(SVector{T,SVector{T,Basic}})
    #jac = zeros(SVector{T,SVector{T,Float64}}) 
    jacDiscrete = zeros(SVector{T,SVector{D,Basic}})
    # inpuVarArr=[]
    u = [symbols("q$i") for i in 1:T] #number of symbols for cont vars
    d = [symbols("d$i") for i in 1:D] #number of symbols for cont vars
    #println("vector of symbols= ",u)
    for i = 3:2+T   # warning the number of diff equ might be less than T 
        jacArr = []
        jacDiscArr = []
        rhs = ((odeExpr.args[i].args[2]))
        # println("rhs= ",(rhs))
        basi = convert(Basic, rhs)#:(2u1+3u2))
        #extract jaco components
        for j = 1:T            # the number of diffequ always coincides with the number of continuous vars?
            coef = diff(basi, u[j])
            #println(typeof(coef))
            push!(jacArr, coef)
        end
        #display(jacArr)
        #extract discreteJaco components
        for j = 1:D            # problem here since we do not know the number of discrete vars
            coef = diff(basi, d[j])
            push!(jacDiscArr, coef)
        end
        #extract inpu vars
        #= for j=1:T
            basi=subs(basi,u[j]=>0)          
        end
        for j=1:D  # this could ve been done in the previous 1:D loop
            basi=subs(basi,d[j]=>0)
        end =#
        jac = setindex(jac, jacArr, i - 2)
        jacDiscrete = setindex(jacDiscrete, jacDiscArr, i - 2)
        # push!(inpuVarArr,basi) 
    end
    #inputVars=SVector{T,Float64}(inpuVarArr)






    Z = M - 2 - T                         #number of if statement ===number of ZC
    ZC_jac = zeros(SVector{Z,SVector{T,Basic}})
    ZC_jacDiscrete = zeros(SVector{Z,SVector{D,Basic}})

    #ZC_inpuVarArr=[]            # temp helper arr
    evsArr = []            # temp helper arr to gather event stucts
    u = [symbols("q$i") for i in 1:T] #number of symbols for cont vars
    zcu = [symbols("x$i") for i in 1:T] #number of symbols for cont vars
    d = [symbols("d$i") for i in 1:D] #number of symbols for disc vars

    for i = 3+T:M   # loop that parses all if statments

        #----------------------------------------------------
        #                       Zero crossing part           #
        #-----------------------------------------------------
        jacArr = []
        jacDiscArr = []
        zc = odeExpr.args[i].args[1].args[2]
        #dump(zc)
        basi = convert(Basic, zc)#:(2u1+3u2))
        #----------------------extract jaco components
        for j = 1:T            # the number of odeProblem always coincides with the number of continuous vars
            coef = diff(basi, zcu[j])
            push!(jacArr, coef)
        end
        ZC_jac = setindex(ZC_jac, jacArr, i - 2 - T)
        #------------------------extract discreteJaco components
        for j = 1:D            # 
            coef = diff(basi, d[j])
            push!(jacDiscArr, coef)
        end
        ZC_jacDiscrete = setindex(ZC_jacDiscrete, jacDiscArr, i - 2 - T)
        #=  #--------------------------extract inpu vars
         for j=1:T
             basi=subs(basi,u[j]=>0)          
         end
         for j=1:D  # this could ve been done in the previous 1:D loop
             basi=subs(basi,d[j]=>0)
         end
         push!(ZC_inpuVarArr,basi) =#
                                                #----------------------------------------------------
                                                #                       events part                 #
                                                #-----------------------------------------------------
        #println("if statement number ",i-2)
        posEvExp = odeExpr.args[i].args[2]  # i corresponds to zc number i; it has 2 events (arg[2]=posEv and arg[3]=NegEv)
        negEvExp = odeExpr.args[i].args[3]
        k = i - 2 - T
        indexPosEv = 2 * k - 1
        indexNegEv = 2 * k
        # now each ev can have many statements 
        #println(length(posEvExp.args))
        #------------------pos Event--------------------#
#=         posEv_disArr = @SVector zeros(D) #  discrete
        posEv_conArr = @SVector zeros(T)  #  continous =#
        posEv_disArr= @SVector fill(NaN, D)
        posEv_conArr= @SVector fill(NaN, T)
        
        for j = 1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
            posEvLHS = posEvExp.args[j].args[1]
            posEvRHS = posEvExp.args[j].args[2]
            #dump(posEvLHS)# symbol at the LHS
            basicLHS = convert(Basic, posEvLHS)

            discVarpositionArray = indexin(basicLHS, d)#basicLHS is a symbol d is a vect of symbols=[d1,d2,d3]    #later try findall(x->x == basicLHS, d)
            #indexin(a, b) Returns a vector containing the highest index in b for each value in a that is a member of b
           # println("discVarpositionArray= ");display(discVarpositionArray);println()
            if !(discVarpositionArray[1] === nothing)
                posEv_disArr = setindex(posEv_disArr, posEvRHS, discVarpositionArray[1])
            else # lhs is not a disc var 
                conVarpositionArray = indexin(basicLHS, u)
                if !(conVarpositionArray[1] === nothing)
                    posEv_conArr = setindex(posEv_conArr, posEvRHS, conVarpositionArray[1])
                else
                    println("LHS is neither a cont nor a discr var!!")
                end
            end
        end

        #------------------neg Event--------------------#
       # negEv_disArr = @SVector zeros(D) #  discrete   # struct DebugUndef x end   #const undef = DebugUndef(NaN)  #missing
        negEv_disArr= @SVector fill(NaN, D)
        #negEv_conArr = @SVector zeros(T)  #  continous
        negEv_conArr=@SVector fill(NaN, T)
        for j = 1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
            negEvLHS = negEvExp.args[j].args[1]
            negEvRHS = negEvExp.args[j].args[2]
            #dump(negEvLHS)# symbol at the LHS
            basicLHS = convert(Basic, negEvLHS)

            discVarpositionArray = indexin(basicLHS, d)

            if !(discVarpositionArray[1] === nothing)
                negEv_disArr = setindex(negEv_disArr, negEvRHS, discVarpositionArray[1])
            else # lhs is not a disc var 
                conVarpositionArray = indexin(basicLHS, u)
                if !(conVarpositionArray[1] === nothing)
                    negEv_conArr = setindex(negEv_conArr, negEvRHS, conVarpositionArray[1])
                else
                    println("LHS is neither a cont nor a discr var!!")
                end
            end
        end

      #println("posEv_disArr= ",posEv_disArr)
        structposEvent = EventDependencyStruct(indexPosEv, posEv_conArr, posEv_disArr)
        push!(evsArr, structposEvent)
        structnegEvent = EventDependencyStruct(indexNegEv, negEv_conArr, negEv_disArr)
        push!(evsArr, structnegEvent)

       #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!there is a problem here:
       #!!!!!!!!!!!!eventDependency: i ignore absent vars by having zero in them, but what about when i actually want to put zero in a var
    end  # end for that parses all if-statments 

    Y = 2 * Z
    eventDependencies = SVector{Y,EventDependencyStruct}(evsArr)  # 2*Z each zc yields 2 events
    #println(eventDependencys)
    # ZCinputVars=SVector{Z,Float64}(ZC_inpuVarArr)
    #println("------------------------ ",jac)
    #based on the type of the problem after a different user input call the appropriate struct
    myodeProblem = NLODEProblem(contVars, discrVars, jac, equs,ZCequs,eventEqus, jacDiscrete, ZC_jac, ZC_jacDiscrete, eventDependencies)
    #myodeProblem=ODEProblem(contVars,discrVars,jac,equs,jacDiscrete,inputVars,ZC_jac,ZC_jacDiscrete,ZCinputVars, eventDependencys)

    myodeProblem






end