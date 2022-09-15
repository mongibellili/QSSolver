
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
    eqs::Expr#Vector{Expr}
    #eqsbang::Vector{Expr}
#=     zceqs::Expr
    eventEqus::Expr =#
    zceqs::Vector{Expr}
    eventEqus::Vector{Expr}
    discreteJacobian::SVector{T,SVector{D,Basic}}
    #inputVars::SVector{T,Float64} 
    ZC_jacobian::SVector{Z,SVector{T,Basic}}
    ZC_jacDiscrete::SVector{Z,SVector{D,Basic}}
    # ZCinputVars::SVector{Z,Float64}   # in case input signal are function of t...use SVector{Z,Function} and create a different struct 
    eventDependencies::SVector{Y,EventDependencyStruct}# 
#=     SD::SVector{T,SVector{T,Int}}
    SZ::SVector{T,SVector{Z,Int}}
    HD::SVector{D,SVector{T,Int}}
    HZ::SVector{D,SVector{Z,Int}} =#
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
    xsymbols=[]############################################################added 9-9-2022
    dsymbols=[]
    cacheSize=12 #default
    casheprovided=false 
    equs=Vector{Expr}()
    num_cache_equs=Vector{Int}()
    num_cache_zcequs=Vector{Int}()
    zcequs=Vector{Expr}()
    eventequs=Vector{Expr}()
    evsArr = []
    postwalk(odeExprs) do x
        if @capture(x, y_ = z_)             # x is expre of type lhs=rhs could ve used head == (:=)
            if y isa Symbol && z isa Integer                    #rhs is int  -> cache provided
                cacheSize=z
                casheprovided=true
            elseif z isa Expr && z.head==:vect # rhs ==vector of state vars initCond or discrete vars
                if y!=:discrete
                    if !initCondreceived
                        initCondreceived=true
                        stateVarName=y
                        du=Symbol(:d,stateVarName)
                        T = length(z.args)
                        contVars = SVector{T,Float64}(z.args) 
                        usymbols = [symbols("q$i") for i in 1:T] # symbols for cont vars
                        xsymbols = [symbols("x$i") for i in 1:T] # symbols for cont vars
                    else 
                        error("initial conditions already defined! for discrete variables, use the identifier 'discrete' for discrete variables")
                    end
                else #y==:discrete #later forbid redefine discrete like cont
                    D = length(z.args)
                    discrVars = Vector{Float64}(z.args)#discrVars = MVector{D,Float64}(z.args)
                     dsymbols = [symbols("d$i") for i in 1:D] # symbols for disc vars
                end
            elseif  du==y.args[1] && ( (z isa Expr && (z.head==:call || z.head==:ref)) || z isa Number)#z is rhs of equations and is an expre:
                #= @show du
                @show z =#
                if z isa Number # rhs of equ =number  not implemented
                    push!(jac, zeros(T))
                    push!(jacDiscrete, zeros(D)) 
                    push!(equs,:($((twoInOneSimplecase(:($(z))))))) 
                    push!(num_cache_equs,1)
                elseif z.head==:ref
                    z=changeVarNames_to_q_d(z,stateVarName)
                    extractJac_from_equs(z,T,D,usymbols,dsymbols,jac,jacDiscrete)
                    push!(equs,:($((twoInOneSimplecase(:($(z))))))) 
                    push!(num_cache_equs,1)
                else

                #cases of z num or ref or expr should be done here and maybe another func like twoinone is implemented
              #=   @show z
                @show :($((twoInOne(:($(z),$(cacheSize))))))
                @show cacheSize
                @show z =#
               # push!(equs,:($((twoInOne(:($(z),$(cacheSize)))).args[1])))    
               z=changeVarNames_to_q_d(z,stateVarName)
               extractJac_from_equs(z,T,D,usymbols,dsymbols,jac,jacDiscrete)
                push!(num_cache_equs,:($((twoInOne(:($(z),$(cacheSize)))).args[2])))   #number of caches distibuted    
                push!(equs,z)  
                end     
               #=  @show num_cache_equs
                @show z   =#        
            else#end of equations
               # error("expression $x: top level contains only expressions 'A=B' or 'if a b' ")#wait until exclude events
            end#end cases inside @capture
        #@show x
       #after capture A=B (init disc equs) we check for 'if sttment'
        elseif x isa Expr && x.head==:if   #@capture if did not work
            !(x.args[1] isa Expr && x.args[1].head==:call && x.args[1].args[1]==:> && (x.args[1].args[3]==0||x.args[1].args[3]==0.0)) && error("use the format 'if a>0")
              x.args[1].args[2] isa Number && error("LHS of term must be be a variable or an expression!")
              x.args[1].args[2]=changeVarNames_to_q_d(x.args[1].args[2],stateVarName)
              extractJac_from_equs(x.args[1].args[2],T,D,usymbols,dsymbols,ZCjac,ZCjacDiscrete)
            if x.args[1].args[2].head==:ref  
                ifexpr=quote
                            function  $(Symbol(:g_, 2))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                                $((twoInOneSimplecase(:($(x.args[1].args[2])))))
                              end 
                         end    
                  push!(zcequs,ifexpr)     
                ########################push!(zcequs,:($((twoInOneSimplecase(:($(x.args[1].args[2]))))))) 
                push!(num_cache_zcequs,1)
            else


                push!(num_cache_zcequs,:($((twoInOne(:($(x.args[1].args[2]),$(cacheSize)))).args[2])))   #number of caches distibuted 
                ifexpr=quote
                    function  $(Symbol(:g_, 2))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
                        $(x.args[1].args[2])
                    end 
                end    
                push!(zcequs,ifexpr)     
           ####################### push!(zcequs,x.args[1].args[2])  
            end     
            #push!(zcequs,:($((twoInOne(:($(x.args[1].args[2]),$(cacheSize)))).args[1])))  
            ##################################################################################################################
            #                                                      events     
            ##################################################################################################################                  
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
    #in the case i am opting for one giant function for all diff equations, equs filled before with RHS of diffequations


    io = IOBuffer() # i guess this is a sys of diff equ solver so I am not optimizing for case 1 equ ; unlike zcf below which can be 0 1 or more
    write(io, "if j==1  $(equs[1]) ;return nothing")
    for i=2:length(equs)-1
        write(io, " elseif j==$i $(equs[i]) ;return nothing")
    end
    write(io, " else  $(equs[length(equs)]) ;return nothing end ")
#=     write(io, "if j==1  $(equs[1]) ;for k=2:$(((num_cache_equs[1]))) cache[k].coeffs.=0.0 end;return nothing")
    for i=2:length(equs)-1
        write(io, " elseif j==$i $(equs[i]) ;for k=2:$(((num_cache_equs[j]))) cache[k].coeffs.=0.0 end;return nothing")
    end
    write(io, " else  $(equs[length(equs)]) ;for k=2:$(((num_cache_equs[length(equs)]))) cache[k].coeffs.=0.0 end;return nothing end ") =#
    s = String(take!(io))
    close(io)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :f   
    def[:args] = [:(j::Int),:(q::Vector{Taylor0{Float64}}),:(d::Vector{Float64}), :(t::Taylor0{Float64}),:(cache::Vector{Taylor0{Float64}})]
    #def[:args] = [:(j::Int),:(derx::Vector{Taylor0{Float64}}),:(q::Vector{Taylor0{Float64}}),:(d::Vector{Float64}), :(t::Taylor0{Float64}),:(cache::Vector{Taylor0{Float64}})]
    def[:body] = Meta.parse(s)
    #def[:rtype]=:nothing# test if cache1 always holds sol  
    functioncode=combinedef(def)

   #=  io = IOBuffer()
    if length(zcequs)>0
        write(io, "if j==1  $(zcequs[1]) ;for k=2:$(((num_cache_zcequs[1]))) cache[k].coeffs.=0.0 end;return nothing")
        for i=2:length(zcequs)-1
            write(io, " elseif j==$i $(zcequs[i]) for k=2:$(((num_cache_zcequs[j]))) cache[k].coeffs.=0.0 end;;return nothing")
        end
        write(io, " else  $(zcequs[length(zcequs)]) ;for k=2:$(((num_cache_zcequs[length(zcequs)]))) cache[k].coeffs.=0.0 end;return nothing end ")
    #= elseif length(zcequs)==1
        write(io, "$(zcequs[1]);return nothing")  =# #LoadError: syntax: "toplevel" expression not at top level      
    else 
        write(io,"return nothing")
    end
    s = String(take!(io))
    close(io)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :zcf   
    def[:args] = [:(j::Int),:(q::Vector{Taylor0{Float64}}),:(d::Vector{Float64}), :(t::Taylor0{Float64}),:(cache::Vector{Taylor0{Float64}})]
    def[:body] = Meta.parse(s)
    zcffunctioncode=combinedef(def)


    io = IOBuffer()
    if length(eventequs)>0 # events come in pairs 0 2 4 6 ....
        write(io, "if j==1  $(eventequs[1]) ;return nothing")
        for i=2:length(eventequs)-1
            write(io, " elseif j==$i $(eventequs[i]) ;return nothing")
        end
        write(io, " else  $(eventequs[length(eventequs)]) ;return nothing end ")
    else
        write(io,"return nothing")
    end
    s = String(take!(io))
    #@show s
    close(io)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = :eventf
    
    def[:args] = [:(j::Int),:(q::Vector{Taylor0{Float64}}),:(d::Vector{Float64}), :(t::Taylor0{Float64}),:(cache::Vector{Taylor0{Float64}})]
    #def[:args] = [:(j::Int),:(derx::Vector{Taylor0{Float64}}),:(q::Vector{Taylor0{Float64}}),:(d::Vector{Float64}), :(t::Taylor0{Float64}),:(cache::Vector{Taylor0{Float64}})]
    def[:body] = Meta.parse(s)
    #def[:rtype]=:nothing# test if cache1 always holds sol
   
    eventfunctioncode=combinedef(def)
 =#

 myodeProblem = NLODEProblem(cacheSize,contVars, discrVars, staticJac, functioncode,zcequs,eventequs, staticjacDiscrete, staticZC_jacobian, staticZC_jacDiscrete, eventDependencies)
    
   # myodeProblem = NLODEProblem(cacheSize,contVars, discrVars, staticJac, functioncode,zcffunctioncode,eventfunctioncode, staticjacDiscrete, staticZC_jacobian, staticZC_jacDiscrete, eventDependencies)
    #myodeProblem=ODEProblem(contVars,discrVars,jac,equs,jacDiscrete,inputVars,ZC_jac,ZC_jacDiscrete,ZCinputVars, eventDependencys)

    myodeProblem
  




    

end