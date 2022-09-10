using MacroTools: postwalk,isexpr,prewalk

#= myregex = r"u\[\d\]"
display(myregex);println()
myexpr=:(u[1]+2)
display(match(myregex,repr(myexpr))!= nothing) =#
myexpr=:(u[1]+2)
dump(myexpr)
println(myexpr.head)
#newexp=postwalk(x -> isexpr(x, Expr) && x.head==:ref  ? x.args[1] : x, myexpr)
#newexp=postwalk(x -> x isa Expr && x.head==:ref  ? Expr(:block,[:($(x.args[1])),:($(x.args[2]))]...) : x, myexpr)
newexp=postwalk(x -> x isa Expr && x.head==:ref  ? Symbol(:u,(x.args[2])) : x, myexpr)
println(myexpr)
println(newexp)
