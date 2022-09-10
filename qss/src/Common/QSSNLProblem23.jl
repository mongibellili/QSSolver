
#helper struct that holds dependency data of an event
struct EventDependencyStruct{T,D} #<: AbstractEventDependecy
    id::Int
    evCont::SVector{T,Float64}
    evDisc::SVector{D,Float64}
end


struct NLODEProblem{T,D,Z,Y}
    cacheSize::Int
    initConditions::SVector{T,Float64}
    #discreteVars::MVector{D,Float64}
    discreteVars::Vector{Float64}
    # jacobian::SVector{T,SVector{T,Float64}} 
    jacobian::SVector{T,SVector{T,Basic}}
    eqs::Vector{Expr}
    #eqsbang::Vector{Expr}
    zceqs::Vector{Expr}
    eventEqus::Vector{Expr}
    discreteJacobian::SVector{T,SVector{D,Basic}}
    #inputVars::SVector{T,Float64} 
    ZC_jacobian::SVector{Z,SVector{T,Basic}}
    ZC_jacDiscrete::SVector{Z,SVector{D,Basic}}
    # ZCinputVars::SVector{Z,Float64}   # in case input signal are function of t...use SVector{Z,Function} and create a different struct 
    eventDependencies::SVector{Y,EventDependencyStruct}# 
end




macro NLodeProblem(odeExprs)
    Base.remove_linenums!(odeExprs)
   # @show odeExprs
    stateVarName=:u #default
    du=:du #default
    initCondreceived=false #bool throw error if user redefine
    #discreteVarName=:d
    T=0#numberStateVars=
    D=0#numberDiscreteVars=
    Z=0#numberzcFunctions=
    contVars=[]
    discrVars=[]
    jac = Vector{Vector{SymEngine.Basic}}()
    jacDiscrete = Vector{Vector{SymEngine.Basic}}()
    ZCjac = Vector{Vector{SymEngine.Basic}}()
    ZCjacDiscrete = Vector{Vector{SymEngine.Basic}}()
    usymbols=[]
    dsymbols=[]
    cacheSize=12 #default
    casheprovided=false 
    equs=Vector{Expr}()
    zcequs=Vector{Expr}()
    eventequs=Vector{Expr}()
    evsArr = []
    postwalk(odeExprs) do x
        if @capture(x, y_ = z_)             # x is expre of type lhs=rhs could ve used head == (:=)
            if y isa Symbol && z isa Integer                    #rhs is int  -> cache provided...should i force user to write "cache"
                cacheSize=z
                casheprovided=true
            elseif z isa Expr && z.head==:vect # rhs ==vector of state vars initCond or discrete vars
                if y!=:discrete
                    if !initCondreceived
                        initCondreceived=true #forbid repetition
                        stateVarName=y
                        du=Symbol(:d,stateVarName)
                        T = length(z.args)
                        contVars = SVector{T,Float64}(z.args) 
                        usymbols = [symbols("q$i") for i in 1:T] # symbols for cont vars
                    else 
                        error("initial conditions already defined for discrete variables!, use the identifier 'discrete' for discrete variables")
                    end
                else #y==:discrete #later forbid redefine discrete like cont
                    D = length(z.args)
                    discrVars = Vector{Float64}(z.args)#discrVars = MVector{D,Float64}(z.args)
                    dsymbols = [symbols("d$i") for i in 1:D] # symbols for disc vars
                end
            elseif  du==y.args[1] && ( (z isa Expr && (z.head==:call || z.head==:ref)) || z isa Number)#z is rhs of diffequations and is an expre or number:
                #= @show du
                @show z =#
                if z isa Number # rhs of equ =number  not implemented
                    push!(jac, zeros(T))
                    push!(jacDiscrete, zeros(D)) 
                else
                    z=changeVarNames_to_q_d(z,stateVarName)
                    extractJac_from_equs(z,T,D,usymbols,dsymbols,jac,jacDiscrete)
                end
                #@show z
                #myex1=(twoInOne(:($(z),$(cacheSize)))).args[1] #modify rhs to fastTaylorFormat #args1 means only eq returned not cache
                fexpr=quote
                    function  $(Symbol(:f_, 1))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                          #  $myex1
#=                           for i=2:10 # first cache holds result
                            cache[i].coeffs.=0.0
                          end =#
                           $((twoInOne(:($(z),$(cacheSize)))).args[1])# warning for users: equations should be given in order, in integrator derx[i] computed                                                                         #from equation number i. if design of integrator changes then we can remove/relax restriction
                        #cache cleanup can be done here: twoinOne return result and caches used, store result in temp, clean cahces used and return temp
                       
                         
                         #$((twoInOne(:($(z),$(cacheSize)))).args[1])
                    end
                    #= function  $(Symbol(:f_, 1))(res::Taylor0{Float64},q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                        #  $myex1
                          res=$((twoInOne(:($(z),$(cacheSize)))).args[1])# warning for users: equations should be given in order, in integrator derx[i] computed                                                                         #from equation number i. if design of integrator changes then we can remove/relax restriction
                  end =#
                end
                push!(equs,fexpr)              
            
            
            else#end of equations
               # error("expression $x: top level contains only expressions 'A=B' or 'if a b' ")#wait until exclude events
            end#end cases inside @capture
        #@show x
       #after capture A=B (init disc equs) we check for 'if sttment'
        elseif x isa Expr && x.head==:if   #@capture if did not work
            !(x.args[1] isa Expr && x.args[1].head==:call && x.args[1].args[1]==:> && (x.args[1].args[3]==0||x.args[1].args[3]==0.0)) && error("use the format 'if a>0")
              # dump(x)
              
              x.args[1].args[2]=changeVarNames_to_q_d(x.args[1].args[2],stateVarName)
              
              extractJac_from_equs(x.args[1].args[2],T,D,usymbols,dsymbols,ZCjac,ZCjacDiscrete)
              #myex1=(twoInOne(:($(z),$(cacheSize)))).args[1] #modify rhs to fastTaylorFormat #args1 means only eq returned not cache
              ifexpr=quote
                  function  $(Symbol(:g_, 1))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                        #  $myex1
                          $((twoInOne(:($(x.args[1].args[2]),$(cacheSize)))).args[1])#                                                            
                  end
              end
             # dump(fexpr)
              push!(zcequs,ifexpr)   
              #for now each pos or neg event has a function...later i can try one event for zc
              x.args[2]=changeVarNames_to_q_d(x.args[2],stateVarName)
              x.args[3]=changeVarNames_to_q_d(x.args[3],stateVarName)
            #=  dump(x.args[2].args[1].args[2])
             dump(x.args[3].args[1].args[2]) =#
            eventexpr=quote
                function  $(Symbol(:g_, 2))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                      #  $myex1
                      #$(x.args[2].args[1].args[1])=$((twoInOne2(:($(x.args[2].args[1].args[2]),$(cacheSize)))).args[1])#  use x.args[2].args[1] to remove outer block  .later deal with d is different than q                                          
                      $(x.args[2])
                    end
            end
            push!(eventequs,eventexpr)  
            eventexpr2=quote
                function  $(Symbol(:g_, 3))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                      #  $myex1
                     # $(x.args[3].args[1].args[1]) =$((twoInOne2(:($(x.args[3].args[1].args[2]),$(cacheSize)))).args[1])#    
                     $(x.args[3])                                            
                 end
            end
            push!(eventequs,eventexpr2)   
              
            #change A[n] to An.
            posEvExp = postwalk(a -> a isa Expr && a.head == :ref ? Symbol((a.args[1]), (a.args[2])) : a, x.args[2])
            negEvExp = postwalk(a -> a isa Expr && a.head == :ref ? Symbol((a.args[1]), (a.args[2])) : a, x.args[3])
              #= posEvExp = x.args[2]  #  'if' has 2 events (arg[2]=posEv and arg[3]=NegEv)
              negEvExp = x.args[3] =#
              indexPosEv = 2 * length(zcequs) - 1 # store events in order
              indexNegEv = 2 * length(zcequs)   
              #------------------pos Event--------------------#
    #=         posEv_disArr = @SVector zeros(D) #  discrete
            posEv_conArr = @SVector zeros(T)  #  continous =#
            posEv_disArr= @SVector fill(NaN, D)   #right now vect of floats
            posEv_conArr= @SVector fill(NaN, T)
            
            for j = 1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
                posEvLHS = posEvExp.args[j].args[1]
                posEvRHS = posEvExp.args[j].args[2]
                #dump(posEvLHS)# symbol at the LHS
                basicLHS = convert(Basic, posEvLHS)

                discVarpositionArray = indexin(basicLHS, dsymbols)#basicLHS is a symbol dsymbols is a vect of symbols=[d1,d2,d3]    #later try findall(x->x == basicLHS, d)
                #indexin(a, b) Returns a vector containing the highest index in b for each value in a that is a member of b
            # println("discVarpositionArray= ");display(discVarpositionArray);println()
                if !(discVarpositionArray[1] === nothing)
                    posEv_disArr = setindex(posEv_disArr, posEvRHS, discVarpositionArray[1])
                else # lhs is not a disc var 
                    conVarpositionArray = indexin(basicLHS, usymbols)
                    if !(conVarpositionArray[1] === nothing)
                        posEv_conArr = setindex(posEv_conArr, posEvRHS, conVarpositionArray[1])
                    else
                        println("LHS is neither a cont nor a discr var!!") #later throw error of used variable not declared in initcond vector or discrete vector
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

                discVarpositionArray = indexin(basicLHS, dsymbols)

                if !(discVarpositionArray[1] === nothing)
                    negEv_disArr = setindex(negEv_disArr, negEvRHS, discVarpositionArray[1])
                else # lhs is not a disc var 
                    conVarpositionArray = indexin(basicLHS, usymbols)
                    if !(conVarpositionArray[1] === nothing)
                        negEv_conArr = setindex(negEv_conArr, negEvRHS, conVarpositionArray[1])
                    else
                        println("LHS is neither a cont nor a discr var!!")
                    end
                end
            end

            structposEvent = EventDependencyStruct(indexPosEv, posEv_conArr, posEv_disArr) # right now posEv_conArr is vect of floats
            push!(evsArr, structposEvent)
            structnegEvent = EventDependencyStruct(indexNegEv, negEv_conArr, negEv_disArr)
            push!(evsArr, structnegEvent)

        #@show x
        end #end cases inside postwalk
           
        



      return x  #if output of pre/postwalk not used then return nothing? not tested
    end #end parent postwalk 

    Z=length(zcequs)
    Y=2*Z
    #later fix error if jac not populated: dimension mismatch

    # zeros(SVector{M,SVector{N,Int}})
  #=   @show jac
    @show size(jac,1)
    @show size(jac,2) =#
    if size(jac,1)==T && length(jac[1])==T
         staticJac = SVector{T,SVector{T,Basic}}(tuple(jac...))
    else
        println("dimension mismatch jac= ",jac)
    end
#=     @show jacDiscrete
    @show size(jacDiscrete,1)
    @show size(jacDiscrete,2) =#

    if size(jacDiscrete,1)==T && length(jacDiscrete[1])==D
        staticjacDiscrete=SVector{T,SVector{D,Basic}}(tuple(jacDiscrete...))
    else
        println("dimension mismatch jacDiscrete= ",jacDiscrete)
    end
    staticZC_jacobian=SVector{Z,SVector{T,Basic}}(tuple(ZCjac...))
    staticZC_jacDiscrete=SVector{Z,SVector{D,Basic}}(tuple(ZCjacDiscrete...))
    
    if !casheprovided
        println("an integer cache size not provided! default cashe used")
    end
      
    
   
    eventDependencies = SVector{Y,EventDependencyStruct}(evsArr)  # 2*Z each zc yields 2 events

    myodeProblem = NLODEProblem(cacheSize,contVars, discrVars, staticJac, equs,zcequs,eventequs, staticjacDiscrete, staticZC_jacobian, staticZC_jacDiscrete, eventDependencies)
    #myodeProblem=ODEProblem(contVars,discrVars,jac,equs,jacDiscrete,inputVars,ZC_jac,ZC_jacDiscrete,ZCinputVars, eventDependencys)

    myodeProblem
  




    

end