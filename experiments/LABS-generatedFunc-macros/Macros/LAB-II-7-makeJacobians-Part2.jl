
module smallqss
using StaticArrays
using SymEngine
using InteractiveUtils


struct Equ{T,D} 
    initConditions ::SVector{T,Float64}  
    discreteVars::MVector{D,Float64} 
    jacobian::SVector{T,SVector{T,Float64}} 
    discreteJacobian::SVector{T,SVector{D,Float64}} 
    inputVars::SVector{T,Float64} 
end
struct ZC_struct{Z,T,D} 
    ZC_jacobian::SVector{Z,SVector{T,Float64}}  
    ZC_discJacobian::SVector{Z,SVector{D,Float64}} 
    ZC_inputVars::SVector{Z,Float64} 
end
struct EventHandlerStruct{T,D} 
    id::Int
    evCont::SVector{T,Float64} 
    evDisc::SVector{D,Float64}
end

 #---------equations----------------------
 macro equations(schema::Expr)
    Base.remove_linenums!(schema)   
    M=length(schema.args) 
    
    T=length((schema.args[1].args[2]).args)
    D=length((schema.args[2].args[2]).args)

   
    contVars=SVector{T,Float64}((schema.args[1].args[2]).args)
    discrVars=MVector{D,Float64}((schema.args[2].args[2]).args)
    
    jac = zeros(SVector{T,SVector{T,Float64}}) 
    jacDiscrete = zeros(SVector{T,SVector{D,Float64}}) 
    inpuVarArr=[]
    u=[symbols("u$i") for i in 1:T] #number of symbols for cont vars
    d=[symbols("d$i") for i in 1:D] #number of symbols for cont vars
    for i=3:M
        jacArr=[]
        jacDiscArr=[]
        rhs=((schema.args[i].args[2]))
        #println(typeof(rhs))
        basi = convert(Basic, rhs)#:(2u1+3u2))
        #extract jaco components
        for j=1:T            # the number of equations always coincides with the number of continuous vars
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
    eq=smallqss.Equ(contVars,discrVars,jac,jacDiscrete,inputVars)
    
 end
 
#---------algorithm--------------------

 macro alg(schema::Expr)
    Base.remove_linenums!(schema)
    #later add code to check user input and throw meaningfull error msg. for example user enter one discrete vars
    #... and writes d2 in event alg
#=     println("-------------begin d1 = 0.1 + u2   end ----------------------")
    println(schema.args[3].args[2])
    println("-------------d1 = 0.1 + u2----------------------")
    println(schema.args[3].args[2].args[1])
    println("-------------d1-------------------------")
    println(schema.args[3].args[2].args[1].args[1]) =#
    M=length(schema.args)       #total number of args: 2 first args for con and disc;the remaining are the number of zc (if statment)
    T=schema.args[1].args[2]    #number of con vars
    D=schema.args[2].args[2]    #number of disc vars
    Z=M-2                         #number of if statement ===number of ZC
    ZC_jac = zeros(SVector{Z,SVector{T,Float64}}) 
    ZC_jacDiscrete = zeros(SVector{Z,SVector{D,Float64}}) 
    
    ZC_inpuVarArr=[]            # temp helper arr
    evsArr=[]            # temp helper arr to gather event stucts
    u=[symbols("u$i") for i in 1:T] #number of symbols for cont vars
    d=[symbols("d$i") for i in 1:D] #number of symbols for disc vars
    for i=3:M   # loop that parses all if statments

                                    #----------------------------------------------------
                                    #                       Zero crossing part           #
                                    #-----------------------------------------------------
        jacArr=[]
        jacDiscArr=[]
        zc=schema.args[i].args[1].args[2] 
        #dump(zc)
        basi = convert(Basic, zc)#:(2u1+3u2))
        #----------------------extract jaco components
        for j=1:T            # the number of equations always coincides with the number of continuous vars
            coef=diff(basi,u[j])
            push!(jacArr,coef)
        end
        ZC_jac=setindex(ZC_jac,jacArr,i-2)
        #------------------------extract discreteJaco components
        for j=1:D            # problem here since we do not know the number of discrete vars...temp solu i made the user add this info before the if statements
            coef=diff(basi,d[j])
            push!(jacDiscArr,coef)
        end
        ZC_jacDiscrete=setindex(ZC_jacDiscrete,jacDiscArr,i-2)
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
        k=i-2
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
     eventHandlers=SVector{2*Z,EventHandlerStruct}(evsArr)  # 2*Z each zc yields 2 events
    #println(eventHandlers)
    ZCinputVars=SVector{Z,Float64}(ZC_inpuVarArr)
    ZC_data=smallqss.ZC_struct(ZC_jac,ZC_jacDiscrete,ZCinputVars)
  (eventHandlers,ZC_data)
    
 end

function solve(myEq::Equ) 
    jac=myEq.jacobian
    u=myEq.initConditions
    d=myEq.discreteVars
    discJac=myEq.discreteJacobian
    inpVar=myEq.inputVars
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

 myEqua=smallqss.@equations(begin
    u=[1.0,2.0,0.5] 
    d=[1.0,0.5]
    du1=u2+2.0     # du1....2 are expected to be in order....later can fix this
    du2=-u1-u2
    du3=u3+u2+d1
end) 
myAlg=smallqss.@alg(begin
    con=3 #number of cont
    dis=2 #number of discr
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
display((myAlg));println() 
#= display((myEqua));println() 
smallqss.solve(myEqua) =#

 #smallqss.solve(myEqua)
