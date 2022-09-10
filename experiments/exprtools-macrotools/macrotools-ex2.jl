using MacroTools: prewalk, postwalk, @capture

#= res= postwalk(x -> x isa Integer ? x + 1 : x, :(2+3))
@show res #:(3 + 4) 
 res= prewalk(x -> x isa Integer ? x + 1 : x, :(2+3))
@show res #:(3 + 4)   =#
##############################################"""
#=res= postwalk(x -> x isa Integer ? :(a + b) : x, 2)
@show res #:(a + b)
 res= prewalk(x -> x isa Integer ? :(a + b) : x, 2)
@show res #:(a + b) =#
##############################################"
#= res= prewalk(x -> x isa Integer ? :(a + b) : x, :(2+3))
@show res #::((a + b) + (a + b)) 
res= postwalk(x -> x isa Integer ? :(a + b) : x, :(2+3))
@show res #::((a + b) + (a + b))=#
########################################"
#= res= postwalk(x -> x isa Integer ? :(a + 1) : x, :(2+3))
@show res #:((a + 1) + (a + 1))
#res= prewalk(x -> x isa Integer ? :(a + 1) : x, :(2+3))
@show res #    error :infinite loop =#

#= ex = quote
    x = f(y, g(z))
    return h(x)
  end
# we want to insert an extra argument 5 into all functions 
  res= postwalk(x -> @capture(x,ff_(args__)) ? :($ff(5,$(args...)) ) : x, ex)
@show res  =#