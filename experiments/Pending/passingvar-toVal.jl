
function createFun(j::Int)
    #f=()->
   
        rhs=j
         f(::Val{j})  =rhs
#=          k=j+1
         f(::Val{k}) where {k} =rhs+1 =#
           
        
   
    
    return f
end

g=createFun(1)
display(g);println()
display(g(Val(1)))