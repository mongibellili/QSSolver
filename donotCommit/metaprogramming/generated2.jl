

using MacroTools: flatten

using BenchmarkTools
#=mutable struct Tupple2FloatList{N, F<:AbstractFloat}
	data :: NTuple{N, NTuple{2, F}}
	index :: Int
	function Tupple2FloatList{N, F}() where {N, F<:AbstractFloat}
		t2fl = new{N, F}()
		t2fl.index = 0
		t2fl
	end
end
#this did not enahnce performance !
@inline function Base.push!(t2fl::Tupple2FloatList{N, F}, val::NTuple{2, F}) where {N, F<:AbstractFloat}
	i = t2fl.index + 1
	GC.@preserve t2fl unsafe_store!(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(t2fl)), convert(NTuple{2, F}, val), i)
	t2fl.index = i
	t2fl
end
=#
@generated function smallest(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	ex = quote m = arg1 < arg2 ? arg1 : arg2 end
	for i in 1:length(args)
		ex = quote
			$ex #value of last m
			 m = args[$i] < m ? args[$i] : m
            #@inbounds m = args[$i] < m ? args[$i] : m
		end
	end
	ex #|> flatten
end
function smallest2(arg1::Float64, arg2::Float64, args::Float64...) 
	 m = arg1 < arg2 ? arg1 : arg2 
	for i in 1:length(args)
		
			
			 m = args[i] < m ? args[i] : m
           
		
	end
	m
end

#display(smallest(1.0, 2.0, -3.0,0.5)) ;println()
x= Array{Float64}(undef,0)
#y=Array{Float64}[]
r=Vector{Array{Float64}}(undef, 1) # it has 1 Array element
#display(r);println()
r[1]=Array{Float64}[]
#push!(r[1],0)
#@show typeof(r[1])
function test()
    for i=1:40000
      # push!(r[1],-i*1.2)
        push!(x,-i*1.2)
       
    end
    
end

#display(x);println()

test()
#display(x)
#@btime smallest(1.0,2.0,x...)
#display(smallest(1.0,2.0,x...))
@btime smallest2(1.0,2.0,x...)
#display(smallest2(1.0,2.0,x...))