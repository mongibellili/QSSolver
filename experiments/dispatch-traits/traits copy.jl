

try_get_single_argtype(f::Function) = map(try_get_single_argtype, methods(f)) # check each method of function

try_get_single_argtype(mm::Method) = try_get_single_argtype(mm.sig)

try_get_single_argtype(x::Any) = nothing # Fail, return nothing to indicate failure.

# We are only intrested in doing things for the 1 argument signature.
# which is the 2-tuple
try_get_single_argtype(::Type{Tuple{F, T}}) where {F, T} = T



is_nothing(::Nothing) = true
is_nothing(::Any) = false

nonscalar_types = unique(filter(!is_nothing, try_get_single_argtype(Base.iterate)))


#println(typeof(nonscalar_types))
nonscalar_types = filter(nonscalar_types) do T
    !occursin(".", string(T)) #Quick hack to see if it is in a module that is loaded
end

#println(nonscalar_types)
aslist(x)=[x]
#= for item in nonscalar_types
    aslist(x::typeof(item))=x
end =#


n= length(nonscalar_types)

#@generated 
#= function popfunc(i::Int)
    
    :(aslist(x::nonscalar_types[$i])=x)
end =#
for i=1:n
    eval(:(aslist(x::nonscalar_types[$i])=x))
   # println(nonscalar_types)
   # aslist(x::nonscalar_types[2])=x
end
#println(popfunc(nonscalar_types))
#= println(aslist([1 2 3]))
println(aslist([1,2,3]))
println(aslist(1))
println(aslist("1"))
 =#
dic=Dict(1=>"one",2=>"two")
println(typeof(dic))
println(aslist(dic))