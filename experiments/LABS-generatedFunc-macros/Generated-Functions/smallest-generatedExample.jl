

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

#= display(smallest(1.0, 2.0, -3.0,0.5))
@btime smallest(1.0, 2.0, -3.0,0.5) =#
display(@code_warntype smallest(1.0, 2.0, -3.0))


function smallest2(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	 m = arg1 < arg2 ? arg1 : arg2 
	for i in 1:length(args)
             m = args[i] < m ? args[i] : m		
	end
	m
end

#= display(smallest2(1.0, 2.0, -3.0,0.5))
@btime smallest2(1.0, 2.0, -3.0,0.5) =#
println("-----------------------------------------------------------------")
display(@code_warntype smallest2(1.0, 2.0, -3.0))