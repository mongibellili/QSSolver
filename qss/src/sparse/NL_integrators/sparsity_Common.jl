#a is a dense jacobian
@inline function getA(::Val{false},cache::MVector{1,Int},a::Vector{Vector{Float64}},i::Int,j::Int,f::Function)# val(1) means prob of type 1 that uses dense jacobians
     a[i][j]
end
#a is a sparse jacobian
@inline function getA(::Val{true},cache::MVector{1,Int},a::Vector{Vector{Float64}},i::Int,j::Int,f::Function)# val(2) means prob of type 1 that uses sparse jacobians through a map function
    f(cache,i,j)
    γ=cache[1]::Int
    if γ!=0
        a[i][γ]::Float64
    else
        0.0
    end
end
@inline function setA(::Val{false},cache::MVector{1,Int},a::Vector{Vector{Float64}},i::Int,j::Int,f::Function,value::Float64)
    a[i][j]=value
    return nothing
end
#= @inline function setotherA(::Val{1},cache::MVector{1,Int},a::Vector{Vector{Float64}},i::Int,j::Int,f::Function,value::Float64)
    a[i][j]=value
    return nothing
end =#
@inline function setA(::Val{true},cache::MVector{1,Int},a::Vector{Vector{Float64}},i::Int,j::Int,f::Function,value::Float64)
    
    f(cache,i,j)
    γ=cache[1]::Int
   #=  if !(j==i || j==i-1 || j==i+1)
    #if !(j in (i,i-1,i+1))
        println("illegal attempt to set a[$i][$j] where j away")
    end =#
    if γ!=0  # the sparse jac structure is fixed upfront, so if for some reasons the code wants to write to an illegal index, it will be ignored...hence in the case equations can change at runtime, user can fill sparse jac with zero upfront (structure aware of future possible fill)
        a[i][γ]=value
    end
    return nothing
end

#= @inline function setotherA(::Val{2},cache::MVector{1,Int},a::Vector{Vector{Float64}},i::Int,j::Int,f::Function,value::Float64)
    
    f(cache,i,j)
    γ=cache[1]::Int
    if !(j==i || j==i-1 || j==i+1)
        println("inside setother: illegal attempt to set a[$i][$j] where j away")
    end
    if γ!=0  # the sparse jac structure is fixed upfront, so if for some reasons the code wants to write to an illegal index, it will be ignored...hence in the case equations can change at runtime, user can fill sparse jac with zero upfront (structure aware of future possible fill)
        a[i][γ]=value
    end
    return nothing
end =#

#= function assignXPrevStepVals(::Val{O},prevStepVal::MVector{T,MVector{O1,Float64}},x::Vector{Taylor0},i::Int)where{T,O,O1}
    #= for m=1:O
        prevStepVal[i][m]=x[i].coeffs[m]
    end =#
    prevStepVal[i] .= x[i].coeffs
end =#


#= function assignXPrevStepVals(::Val{O},prevStepVal::Vector{Float64} ,x::Vector{Taylor0},i::Int)where{O}
    
        prevStepVal[i]=x[i][0]
   
    #prevStepVal[i] .= x[i].coeffs
end =#

#= 
 function getPrevStepVal(prevStepVal::MVector{T,MVector{O1,Float64}},index::Int)where{T,O1}
              prevStepVal[index][1]
 end
 function getPrevStepVal(prevStepVal::Vector{Float64},index::Int)
    prevStepVal[index]
end =#