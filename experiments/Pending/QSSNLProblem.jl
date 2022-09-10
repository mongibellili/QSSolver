
#helper struct that holds dependency data of an event
struct EventDependencyStruct{T,D} #<: AbstractEventDependecy
    id::Int
    evCont::SVector{T,Float64}
    evDisc::SVector{D,Float64}
end


struct NLODEProblem{T,D,Z,Y}
    initConditions::SVector{T,Float64}
    #discreteVars::MVector{D,Float64}
    discreteVars::Vector{Float64}
    # jacobian::SVector{T,SVector{T,Float64}} 
    jacobian::SVector{T,SVector{T,Basic}}
    eqs::Vector{Expr}
    eqsbang::Vector{Expr}
    zceqs::Vector{Expr}
    eventEqus::Vector{Expr}
    discreteJacobian::SVector{T,SVector{D,Basic}}
    #inputVars::SVector{T,Float64} 
    ZC_jacobian::SVector{Z,SVector{T,Basic}}
    ZC_jacDiscrete::SVector{Z,SVector{D,Basic}}
    # ZCinputVars::SVector{Z,Float64}   # in case input signal are function of t...use SVector{Z,Function} and create a different struct 
    eventDependencies::SVector{Y,EventDependencyStruct}# i think it can better if use matrix instead of EventDependencyStruct
end


struct stub
    x::Int
end
function changeVarNames_to_q_d(ex::Expr,stateVarName::Symbol)
    newEx=postwalk(ex) do a  # change rhs equ to contain q[] and d[] instead of user symbols
        if a isa Expr && a.head == :ref # a is expre of type A[n]  ie a=A[n]
            a.args[1]!=stateVarName && a.args[1]!=:discrete && error("symbol $(a.args[1]) undefined")
            if a.args[1]==stateVarName 
               a.args[1]=:q
            else
                a.args[1]=:d
            end
        end
        return a
    end
    newEx
end
function extractJac_from_equs(ex::Expr,T::Int,D::Int,usymbols::Vector{SymEngine.Basic},dsymbols::Vector{SymEngine.Basic},jac::Vector{Vector{SymEngine.Basic}},jacDiscrete::Vector{Vector{SymEngine.Basic}})
    m=postwalk(ex) do a   #after equs constructed, eliminate ref ; use new expr m since we need z for equs below
        if a isa Expr && a.head == :ref # a is expre of type A[n]  ie a=A[n]
        a=Symbol((a.args[1]), (a.args[2]))  #a become An #needed for differentiation ...jacobians....
        end
        return a
    end
    #@show m  #m = :((q2 - q1) + 3.2)
    jacArr = []
    jacDiscArr = []
    basi = convert(Basic, m)
    #extract jaco components
    for j = 1:T            
        coef = diff(basi, usymbols[j])
        push!(jacArr, coef)
    end
    for j = 1:D            
        coef = diff(basi, dsymbols[j])
        push!(jacDiscArr, coef)
    end
    push!(jac, jacArr)
    push!(jacDiscrete, jacDiscArr)             
end

macro NLodeProblem(odeExprs)
    Base.remove_linenums!(odeExprs)
    #default vars
   # @show odeExprs
    stateVarName=:u
    initCondreceived=false #bool throw error if user redefine
    #discreteVarName=:d
    T=0#numberStateVars=
    D=0#numberDiscreteVars=
    contVars=[]
    discrVars=[]
    jac = Vector{Vector{SymEngine.Basic}}()
    jacDiscrete = Vector{Vector{SymEngine.Basic}}()
    usymbols=[]
    dsymbols=[]
    casheSize=12 #default
    casheprovided=false
   
    equs=Vector{Expr}()
    postwalk(odeExprs) do x
       if @capture(x, y_ = z_)             # x is expre of type lhs=rhs
            if z isa Integer                    #rhs is int  -> cache provided
                casheSize=z
                casheprovided=true
            elseif z isa Expr && z.head==:vect # rhs ==vector of state vars initCond or discrete vars
                if y!=:discrete
                    if !initCondreceived
                        initCondreceived=true
                        stateVarName=y
                        T = length(z.args)
                        contVars = SVector{T,Float64}(z.args) 
                        #jac = zeros(SVector{T,SVector{T,Basic}})
                        usymbols = [symbols("q$i") for i in 1:T] # symbols for cont vars
                    else 
                        error("initial conditions already defined! for discrete variables, use the identifier 'discrete' ")
                    end
                else #y==:discrete #later forbid redefine discrete like cont
                    D = length(z.args)
                    discrVars = Vector{Float64}(z.args)
                    #jacDiscrete = zeros(SVector{T,SVector{D,Basic}})
                     #discrVars = MVector{D,Float64}(z.args)
                     dsymbols = [symbols("d$i") for i in 1:D] # symbols for disc vars
                end
            elseif z isa Expr && z.head==:call #z is rhs of equations
                z=changeVarNames_to_q_d(z,stateVarName)
                #@show z
                #= m=postwalk(z) do a   #after equs constructed, eliminate ref ; use new expr m since we need z for equs below
                    if a isa Expr && a.head == :ref # a is expre of type A[n]  ie a=A[n]
                       a=Symbol((a.args[1]), (a.args[2]))  #a become An #needed for differentiation ...jacobians....
                    end
                    return a
                end
                #@show m  #m = :((q2 - q1) + 3.2)
                jacArr = []
                jacDiscArr = []
                basi = convert(Basic, m)
                #extract jaco components
                for j = 1:T            
                    coef = diff(basi, usymbols[j])
                    push!(jacArr, coef)
                end
                for j = 1:D            
                    coef = diff(basi, dsymbols[j])
                    push!(jacDiscArr, coef)
                end
                push!(jac, jacArr)
                push!(jacDiscrete, jacDiscArr)              =#  
                extractJac_from_equs(z,T,D,usymbols,dsymbols,jac,jacDiscrete)


                #myex1=(twoInOne(:($(z),$(casheSize)))).args[1] #modify rhs to fastTaylorFormat #args1 means only eq returned not cache
                fexpr=quote
                    function  $(Symbol(:f_, 1))(q::Vector{Taylor1{Float64}},d::Vector{Float64}, t::Taylor1{Float64},cache::Vector{Taylor1{Float64}})
                          #  $myex1
                            $((twoInOne(:($(z),$(casheSize)))).args[1])
                    end
                end
                push!(equs,fexpr)              
            #end of equations
            else
                println("expression $x: $y=$z")
                error("problem contains expression 'A=B' that does not meet the required problem specs! ")
            end#end cases inside @capture
        
       #after capture A=B (init disc equs) we check for 'if sttment'
        elseif x isa Expr && x.head==:if 
            !(x.args[1] isa Expr && x.args[1].head==:call && x.args[1].args[1]==:> && (x.args[1].args[3]==0||x.args[1].args[3]==0.0)) && error("use the format 'if a>0")
            
        end #end cases inside postwalk
           
        



      return x  #if output of pre/postwalk not used then return nothing? not tested
    end #end parent postwalk 

    staticJac = SVector{T,SVector{T,Basic}}(tuple(jac...))
    staticjacDiscrete=SVector{T,SVector{D,Basic}}(tuple(jacDiscrete...))
  #=  @show newEX
    @show odeExprs =#

  #=   @show discrVars
    @show contVars

    @show stateVarName
    @show T
    @show D
    @show contVars
    @show discrVars =#
   
   @show staticJac
    @show staticjacDiscrete 
   #=  @show usymbols
    @show dsymbols
    @show equs =#
    if !casheprovided
        println("an integer cache size not provided! default cashe used")
    end
      
    
   

  




    


mystub=stub(2) 
mystub
end