

using MacroTools: flatten
using InteractiveUtils
using BenchmarkTools

@generated function smallest(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	ex = quote m = arg1 < arg2 ? arg1 : arg2 end
	for i in 1:length(args)
		ex = quote
			$ex #value of last m
			#@inbounds m = args[$i] < m ? args[$i] : m
             m = args[$i] < m ? args[$i] : m
		end
	end
	ex #|> flatten
end

#display(smallest(1.0, 2.0, -3.0,0.5))#-3.0
#@btime smallest(1.0, 2.0, -3.0,0.5) #0.014 ns (0 allocations: 0 bytes)
#display(@code_warntype smallest(1.0, 2.0, -3.0))


function smallest2(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	 m = arg1 < arg2 ? arg1 : arg2 
	for i in 1:length(args)
             m = args[i] < m ? args[i] : m		
	end
	m
end

display(smallest2(1.0, 2.0, -3.0,0.5))#-3.0
@btime smallest2(1.0, 2.0, -3.0,0.5)#1.646 ns (0 allocations: 0 bytes)
println("-----------------------------------------------------------------")
#display(@code_warntype smallest2(1.0, 2.0, -3.0))
#=the body of the generated function was only executed once here, for the specific set of argument types, and the result was cached
	After that, the expression returned from the generated function on the first invocation is re-used. if all the information we need 
	 is embedded in the type information of the arguments, we can utilize generated functions to move the iteration to compile-time. then we have a simpler expression
	 than the original code that is returned evey time, and evaluated at run time.
	warning: if a code not in exprssion it would not run after the first invocation.
	=#