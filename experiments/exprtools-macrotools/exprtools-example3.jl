#using ExprTools
using MacroTools
macro log_trace(expr)
    def=splitdef(expr)
   # display(def);println()
   args, kwargs, arg_dict = get_args(func_def)
   display(args);println()
   display(kwargs);println()
   display(arg_dict);println()

   # push!(def,:name => :g)
    name=Meta.quot(get(def,:name,Symbol("<anon>")))
    def[:body]=quote
        #println("entering ", $name, $(args_tuple_expr(def)))# any code that refers x.. necessitates existence of combinedef !!! looks like x is not defined unless combinedef is written
       # $(def[:body])
       println("hello")# this does not require combinedef
    end
   
   # def[:name] = :g
  
   #display(def)
    combinedef(def)

end
a=@log_trace function foo(x,y)
    return y*x
end
#display(foo(2,3))
#display((methods(g)))
#=  @log_trace bar(x)=3*x
 #display((methods(g)))#2 methods for generic function "g"
@log_trace x->4*x
#display(qux)   =#
#= foo(1)
bar(2)
qux(3) =#