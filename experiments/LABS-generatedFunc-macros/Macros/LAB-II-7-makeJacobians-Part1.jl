
module smallqss
using StaticArrays
using SymEngine
using InteractiveUtils

struct VarsStruct{T,D} 
    initConditions ::SVector{T,Float64}  
    discreteVars::SVector{D,Float64} 
end
struct Equ{T} 
    jacobian::SVector{T,SVector{T,Float64}} 
end
struct Alg{T} 
    ZC_jacobian::SVector{T,Float64}   
end
#---------vars------------------------
macro varis(schema::Expr)
    Base.remove_linenums!(schema)
    #cont=:($schema).args[1].args[2]
    (quote  
    N=length($(schema.args[1].args[2])) 
    D=length($(schema.args[2].args[2])) 
    u=SVector{N,Float64}($(schema.args[1].args[2]))
    d=SVector{D,Float64}($(schema.args[2].args[2]))
    var=smallqss.VarsStruct(u,d)
    end)

end

 #---------equations----------------------
 macro equations(schema::Expr)
    Base.remove_linenums!(schema)   
    N=length(schema.args) 
    jac = zeros(SVector{N,SVector{N,Float64}}) 
    #println(N)
    u=[symbols("u$i") for i in 1:N]
    for i=1:N
        jacArr=[]
        rhs=((schema.args[i].args[2]))
        #println(typeof(rhs))
        for j=1:N
            basi = convert(Basic, rhs)#:(2u1+3u2))
            coef=diff(basi,u[j])
            push!(jacArr,coef)
        end
        jac=setindex(jac,jacArr,i)
    end  
    eq=smallqss.Equ(jac)
 end
 #=
#---------algorithm--------------------
 function alg_to_Jacs(schema)
   # println(schema)   
         esc(quote                        
                end)
 end
 macro alg(schema::Expr)
     Base.remove_linenums!(schema)
     rhs=:($schema).args[1].args[2]
    # println(rhs)
     alg_to_Jacs(rhs) 
 end =#


function testVars(myVars::VarsStruct)
    display(myVars.initConditions[2]+3.2);println() 
end


function solve(myEq::Equ,myVars::VarsStruct) 
    jac=myEq.jacobian
    u=myVars.initConditions
       # for index=1:length(u)
       index=3
            computeDerivative(index,jac,u)
       # end
end 
function computeDerivative(index::Int,jacobian::SVector{T,SVector{T,Float64}},u ::SVector{T,Float64}  ) where {T} 
       der= 0.0
       for j = 1:T
       der += jacobian[index][j] * u[j] 
       end
       display(der);println()
end 

end #end module
#--------------------------user space--------------------
using StaticArrays
using BenchmarkTools
myVars=smallqss.@varis(begin
    u=[1.0,2.0,0.5] 
    d=[0.0]
end) 
 myEqua=smallqss.@equations(begin
    du[1]=u2 
    du[2]=-u1-u2
    du[3]=u3+u2
end) 
#=myAlg=smallqss.@alg(begin
    if u[1] >0 
        d[1]=0.0
    else
        d[1]=1.0
    end

end) =#

display((myEqua));println() 
#= display(dump(myVars));println() 
display(typeof(myVars));println()  =#
#smallqss.testVars(myVars)
smallqss.solve(myEqua,myVars)
#display(u)
#@btime smallqss.solve(u)
