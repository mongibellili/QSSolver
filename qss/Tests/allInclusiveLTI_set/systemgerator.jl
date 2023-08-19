using LinearAlgebra
struct LTIProblem
    name::String
    x0::Vector{Float64}#initconds
    input::Vector{Float64}
    jac::Array{Float64, 2}
end
struct LTISolution
    name::String
    eignVal::Vector{Float64}
    eignVec::Vector{Float64}
    solPart::Vector{Float64}
    coefs::Vector{Float64}
end
function generateSol(ltiprob::LTIProblem)
    A=ltiprob.jac;b=ltiprob.input;x10=ltiprob.x0[1];x20=ltiprob.x0[2]
    xp=-inv(A)*b 
    V=eigvecs(A)
    λ=eigvals(A)
   # @show λ[1],λ[2]
    V1=V[1]/V[2] 
    V2=V[3]/V[4] 
  #  @show V1,V2
    c2=(x10-xp[1]-(x20-xp[2])*V1)/(V2-V1)
    c1=x20-xp[2]-c2
   #=  @show xp[1],xp[2]
    @show c1,c2 =#
    LTISolution(ltiprob.name,[λ[1],λ[2]],[V1,V2],[xp[1],xp[2]],[c1,c2])
end
function consructJacobiansSystemC(coeffs:: Vector{Float64},coeffsSpec:: Vector{Float64})
     allJacsysC=Vector{ Array{Float64, 2}}()
     for a11 in coeffsSpec
        if a11<0
            for a22 in coeffs
                if a22<0
                    for a12 in coeffs
                        for a21 in coeffs
                            if a12*a21>0#same sign
                                if a11*a22>a12*a21
                                    A=[a11 a12;a21 a22]
                                    push!(allJacsysC,A)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
   #  push!(allJacsysC,A)
    # display(allJacsysC)
    return allJacsysC
end
function writeProblemstoFile()
    coefs=[-20.0,-1.0,1.0,20.0]
    coefsSpec=[-21.0,-1.1,1.0,20.0]
    #coefs=[-1.0,1.0];coefsSpec=[-1.1,1.0]
    jacs=consructJacobiansSystemC(coefs,coefsSpec)
   # allCprobs=Vector{LTIProblem}(); allCsols=Vector{LTISolution}()
    initConds=[-1.0,-2.0]
    us=1;ul=20
    inputs=[[us,us],[-us,-us],[-us,us],[us,ul],[-us,-ul],[-us,ul],[us,-ul],[ul,ul],[-ul,-ul],[-ul,ul]]
   # inputs=[[us,us]]
    counter=0
    path="./Tests/typeCCC.jl" #default path
    vectprStr="["  # helper to manually written vector later in maintest
    for jac in jacs
        for input in inputs
            counter+=1
            lti_prob=LTIProblem("C$counter",initConds,input,jac)
           # push!(allCprobs,lti_prob) 
            anaSol=generateSol(lti_prob)
           # push!(allCsols,anaSol)
          
           vectprStr*="C$(counter),"
    
            ss="\n function $(lti_prob.name)() \n"
            ss*=" odeprob = @NLodeProblem begin
                name=($(lti_prob.name),)
                u = [$(lti_prob.x0[1]), $(lti_prob.x0[2])]
                du[1] = $(lti_prob.jac[1,1])*u[1]+$(lti_prob.jac[1,2])*u[2]+$(lti_prob.input[1])
                du[2] =$(lti_prob.jac[2,1])*u[1]+$(lti_prob.jac[2,2])*u[2]+$(lti_prob.input[2])
            end  "
            ss*="\n x1(t)=$(anaSol.coefs[1])*$(anaSol.eignVec[1])*exp($(anaSol.eignVal[1])*t)+$(anaSol.coefs[2])*$(anaSol.eignVec[2])*exp($(anaSol.eignVal[2])*t)+$(anaSol.solPart[1])"
            ss*="\n x2(t)=$(anaSol.coefs[1])*exp($(anaSol.eignVal[1])*t)+$(anaSol.coefs[2])*exp($(anaSol.eignVal[2])*t)+$(anaSol.solPart[2])"
            ss*="\n return (odeprob,x1,x2) \n end "

            open(path, "a") do io    
                println(io,ss) 

            end

        end
    end
    vectprStr*="]"
    open(path, "a") do io    
        println(io,vectprStr) 

    end
end
writeProblemstoFile()