using MacroTools
println("----------------@capture-----------------")

ex = quote
    struct Foo1
      x::Int
      y
    end
  end
  #dump(ex.args[2])
 # println(isexpr(ex.args[2], :struct)) #is head of exprssion in (:struct,...)  #true
#=  @capture(ex, struct T_ fields__ end) #return bool and mute vars store exprs
#@show T#:Foo
#@show fields# [:(x::Int), :y]

ex1=:( bar(x::N) where {N<:Number} = 3 )
@capture(ex1, f_(xs__) where {T_} = body_) ||
  error("expected a function with a single type parameter") =#
#@show f,xs,T,body  #(:bar, Any[:(x::N)], :(N <: Number), quote #= /home/unknown/QS_Solver/experiments/exprtools-macrotools/macrotools-ex1.jl:18 =#  3   end)



 #capture(:(foo("$(a)")), foo(x_))#error foo not defined
# capture(:(foo("$(a)")), foo(x_String))##error foo not defined
#=   @capture(:(foo("$(a)ghgfh")), foo(x_String_string)) #this matches expre by the type
  @show x #:("$(a)")
  @capture(:(foo("4")), foo(x_String_string))#nothing for 4 :(4) :("4")
  @show x#"4" =#
  @capture(:(foo(5)), foo(x_Int_f)) #this matches expre by the type  #require write something after type???!!!!
  #@capture(:(foo(5.0)), foo(x_Int_f)) #nothing
  @show x #:("$(a)")
# f=Symbol(:t)
 # @capture(:(foo(f)), foo(x_Symbol)) #the symbol doesnot require to add something after it

 
  #@show x #:("$(a)")