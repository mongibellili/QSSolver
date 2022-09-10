
module smallqss
using StaticArrays
using SymEngine
using InteractiveUtils


#abstract type  AbstractEventHandler end
struct EventHandlerStruct{T,D} #<: AbstractEventHandler
    id::Int
    evCont::SVector{T,Float64} 
    evDisc::SVector{D,Float64}
end



struct ODEProblem{T,D,Z,Y} 
    initConditions ::SVector{T,Float64}  
    discreteVars::MVector{D,Float64} 
    jacobian::SVector{T,SVector{T,Float64}} 
    discreteJacobian::SVector{T,SVector{D,Float64}} 
    inputVars::SVector{T,Float64} 
    ZC_jacobian::SVector{Z,SVector{T,Float64}}  
    ZC_jacDiscrete::SVector{Z,SVector{D,Float64}} 
    ZCinputVars::SVector{Z,Float64} 
    eventHandlers::SVector{Y,EventHandlerStruct}
end








 #---------odeProblem----------------------
 macro odeProblem(schema::Expr)
    Base.remove_linenums!(schema)   
    M=length(schema.args)   # total number of args... 1 for initivars+1 for discretevars + T for num_diff_eqs + M-2-T for ZC if_statments
    
    T=length((schema.args[1].args[2]).args)  # number of cont vars.....correponds to the number of diff equ 
    D=length((schema.args[2].args[2]).args)   # number of discrete vars

                                    #--------------------------------------------------------
                                    #               intial vals and diff equa part           #
                                    #---------------------------------------------------------
   
    contVars=SVector{T,Float64}((schema.args[1].args[2]).args)
    discrVars=MVector{D,Float64}((schema.args[2].args[2]).args)
    
    jac = zeros(SVector{T,SVector{T,Float64}}) 
    jacDiscrete = zeros(SVector{T,SVector{D,Float64}}) 
    inpuVarArr=[]
    u=[symbols("u$i") for i in 1:T] #number of symbols for cont vars
    d=[symbols("d$i") for i in 1:D] #number of symbols for cont vars
    for i=3:2+T
        jacArr=[]
        jacDiscArr=[]
        rhs=((schema.args[i].args[2]))
        #println(typeof(rhs))
        basi = convert(Basic, rhs)#:(2u1+3u2))
        #extract jaco components
        for j=1:T            # the number of odeProblem always coincides with the number of continuous vars
            coef=diff(basi,u[j])
            push!(jacArr,coef)
        end
        #extract discreteJaco components
        for j=1:D            # problem here since we do not know the number of discrete vars
            coef=diff(basi,d[j])
            push!(jacDiscArr,coef)
        end
        #extract inpu vars
        for j=1:T
            basi=subs(basi,u[j]=>0)          
        end
        for j=1:D  # this could ve been done in the previous 1:D loop
            basi=subs(basi,d[j]=>0)
        end
        jac=setindex(jac,jacArr,i-2)
        jacDiscrete=setindex(jacDiscrete,jacDiscArr,i-2)
        push!(inpuVarArr,basi)
    end  
    inputVars=SVector{T,Float64}(inpuVarArr)
   





    Z=M-2-T                         #number of if statement ===number of ZC
    ZC_jac = zeros(SVector{Z,SVector{T,Float64}}) 
    ZC_jacDiscrete = zeros(SVector{Z,SVector{D,Float64}}) 
    
    ZC_inpuVarArr=[]            # temp helper arr
    evsArr=[]            # temp helper arr to gather event stucts
    u=[symbols("u$i") for i in 1:T] #number of symbols for cont vars
    d=[symbols("d$i") for i in 1:D] #number of symbols for disc vars
    for i=3+T:M   # loop that parses all if statments

                                    #----------------------------------------------------
                                    #                       Zero crossing part           #
                                    #-----------------------------------------------------
        jacArr=[]
        jacDiscArr=[]
        zc=schema.args[i].args[1].args[2] 
        #dump(zc)
        basi = convert(Basic, zc)#:(2u1+3u2))
        #----------------------extract jaco components
        for j=1:T            # the number of odeProblem always coincides with the number of continuous vars
            coef=diff(basi,u[j])
            push!(jacArr,coef)
        end
        ZC_jac=setindex(ZC_jac,jacArr,i-2-T)
        #------------------------extract discreteJaco components
        for j=1:D            # 
            coef=diff(basi,d[j])
            push!(jacDiscArr,coef)
        end
        ZC_jacDiscrete=setindex(ZC_jacDiscrete,jacDiscArr,i-2-T)
        #--------------------------extract inpu vars
        for j=1:T
            basi=subs(basi,u[j]=>0)          
        end
        for j=1:D  # this could ve been done in the previous 1:D loop
            basi=subs(basi,d[j]=>0)
        end
        push!(ZC_inpuVarArr,basi)
                                    #----------------------------------------------------
                                    #                       events part                 #
                                    #-----------------------------------------------------
        #println("if statement number ",i-2)
        posEvExp=schema.args[i].args[2]  # i corresponds to zc number i; it has 2 events (arg[2]=posEv and arg[3]=NegEv)
        negEvExp=schema.args[i].args[3]
        k=i-2-T
        indexPosEv=2*k-1
        indexNegEv=2*k
        # now each ev can have many statements 
        #println(length(posEvExp.args))
                         #------------------pos Event--------------------#
        posEv_disArr=@SVector zeros(D) #  discrete
        posEv_conArr=@SVector zeros(T)  #  continous
        for j=1:length(posEvExp.args)  # j coressponds the number of statements under one posEvent
           posEvLHS=posEvExp.args[j].args[1]
           posEvRHS=posEvExp.args[j].args[2]
           #dump(posEvLHS)# symbol at the LHS
           basicLHS = convert(Basic, posEvLHS)
           
           discVarpositionArray = indexin(basicLHS, d)
           
           if !(discVarpositionArray[1] === nothing)
            posEv_disArr=setindex(posEv_disArr,posEvRHS,discVarpositionArray[1])
           else # lhs is not a disc var 
                conVarpositionArray = indexin(basicLHS, u)
                if !( conVarpositionArray[1]=== nothing)
                    posEv_conArr=setindex(posEv_conArr,posEvRHS,conVarpositionArray[1])
                else
                        println("LHS is neither a cont nor a discr var!!")
                end
            end
        end

              #------------------neg Event--------------------#
              negEv_disArr=@SVector zeros(D) #  discrete
              negEv_conArr=@SVector zeros(T)  #  continous
              for j=1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
                 negEvLHS=negEvExp.args[j].args[1]
                 negEvRHS=negEvExp.args[j].args[2]
                 #dump(negEvLHS)# symbol at the LHS
                 basicLHS = convert(Basic, negEvLHS)
                 
                 discVarpositionArray = indexin(basicLHS, d)
                 
                 if !(discVarpositionArray[1] === nothing)
                  negEv_disArr=setindex(negEv_disArr,negEvRHS,discVarpositionArray[1])
                 else # lhs is not a disc var 
                      conVarpositionArray = indexin(basicLHS, u)
                      if !( conVarpositionArray[1]=== nothing)
                          negEv_conArr=setindex(negEv_conArr,negEvRHS,conVarpositionArray[1])
                      else
                              println("LHS is neither a cont nor a discr var!!")
                      end
                  end
              end


       structposEvent=smallqss.EventHandlerStruct(indexPosEv,posEv_conArr,posEv_disArr)
       push!(evsArr,structposEvent)
       structnegEvent=smallqss.EventHandlerStruct(indexNegEv,negEv_conArr,negEv_disArr)
       push!(evsArr,structnegEvent)
     

    end  # end for that parses all if-statments 
     #------------------instantiate the struct
     #eq=smallqss.Equ(contVars,discrVars,jac,jacDiscrete,inputVars)
     Y=2*Z
     eventHandlers=SVector{Y,EventHandlerStruct}(evsArr)  # 2*Z each zc yields 2 events
    #println(eventHandlers)
    ZCinputVars=SVector{Z,Float64}(ZC_inpuVarArr)
    #ZC_data=smallqss.ZC_struct(ZC_jac,ZC_jacDiscrete,ZCinputVars)

    myodeProblem=smallqss.ODEProblem(contVars,discrVars,jac,jacDiscrete,inputVars,ZC_jac,ZC_jacDiscrete,ZCinputVars, eventHandlers)

    myodeProblem
    
    
 end
 

function solve(mypr::ODEProblem) 
    jac=mypr.jacobian
    u=mypr.initConditions
    d=mypr.discreteVars
    discJac=mypr.discreteJacobian
    inpVar=mypr.inputVars
    zc_jac=mypr.ZC_jacobian
    ZC_discJac=mypr.ZC_jacDiscrete
    ZC_input=mypr.ZCinputVars
    evHandlr=mypr.eventHandlers
#=     println(ZC_discJac)
    println(ZC_input)
    println(evHandlr) =#
    #for k=1:100
        for index=1:length(u)
      # index=3
            computeDerivative(index,jac,u,d,discJac,inpVar)

        end
    #end
end 
function computeDerivative(index::Int,jacobian::SVector{T,SVector{T,Float64}},u ::SVector{T,Float64} , d::MVector{D,Float64},  discJac::SVector{T,SVector{D,Float64}} ,inputVars::SVector{T,Float64} ) where {T,D} 
       der= 0.0
       for j = 1:T
         der += jacobian[index][j] * u[j] 
       end
       for j = 1:D
         der += discJac[index][j]*d[j]
       end
       der+=  inputVars[index]
       display(der);println()
end 

end #end module
#--------------------------user space--------------------
using StaticArrays
using BenchmarkTools

 myProb=smallqss.@odeProblem(begin
    u=[1.0,2.0,0.5] 
    d=[1.0,0.5]
    du1=u2+2.0     # du1....du n are expected to be in order....later can fix this
    du2=-u1-u2
    du3=u3+u2+d1
    if u1+0.7 >0 
        d1=0.1
    else
        d1=1.0
    end
    if u2+d1 >0 
        d2=0.33
        u3=2.2
    else
        d2=1.0
    end
end) 


display(myProb);println() 

#= display((myEqua));println() 
smallqss.solve(myEqua) =#

 #smallqss.solve(myProb)
