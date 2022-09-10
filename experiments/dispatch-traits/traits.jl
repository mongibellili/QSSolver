#println(String isa Type{<:String}) #true
#println( Type{<:String} )
#println(Type{<:String} isa Type{String}  ) #false
#println(Type{<:String} isa Type  )#true
#println(DataType <: Type{Int})#false
#println(DataType <: Type{T} where{T})#true
#println(supertype(DataType)===Type{T} where{T})#false
#println(supertype(DataType))#Type{T}
#println(supertype(UnionAll))#Type{T}
#println(supertype(Type))#Any
#println(supertype(Any))#Any
#println(typeof(Type{T}where{T}))#error T not defined
#println(typeof(Type))#UnionAll
#println(typeof(UnionAll))#DataType
#println(typeof(DataType))#DataType
#println(typeof(Any))#DataType






aslist_direct(x) = [x]
#aslist_direct(x::Union{AbstractArray,Tuple}) = x


#= @show aslist_direct([1,2,3])#[1, 2, 3]
@show aslist_direct(1);#[1]
 =#
# but wait do we have to define aslist_direct for every type we know
#yes otherwise string will boxed

abstract type Scalarness end
struct Scalar <: Scalarness end
struct NonScalar <: Scalarness end

#in dispatch (previous file) i used this newly created struct type to do dispatch on in the arguments
#traits does something different: instances of this newly created struct types are used as the RHS for 
# a trait-function: so this function returns our new type

scalarness(::Type) = Scalar() #fall-back, by default everything is scalar
scalarness(::Type{<:AbstractArray}) = NonScalar()
scalarness(::Type{<:Tuple}) = NonScalar()
scalarness(::Type{<:AbstractString}) = NonScalar()

# other types instead of used in arguments they are used in the trait function which is put in the arguments
# now the dispacth will be on the return type of that function

aslist(x::T) where T = aslist(scalarness(T), x)
aslist(::Scalar, x) = [x]
aslist(::NonScalar, x) = x

# this requires us also to define the function scalarness on types we know!!
#but it is still beneficial since if we need another functionality like aslist we
#do not need to redo the whole work

#= 
using InteractiveUtils
println(@edit methods(Base.iterate)) =#
#= sig=methods(Base.iterate).ms[2].sig
println(sig)
 =#



try_get_single_argtype(f::Function) = map(try_get_single_argtype, methods(f)) # check each method of function

try_get_single_argtype(mm::Method) = try_get_single_argtype(mm.sig)

try_get_single_argtype(x::Any) = nothing # Fail, return nothing to indicate failure.

# We are only intrested in doing things for the 1 argument signature.
# which is the 2-tuple
try_get_single_argtype(::Type{Tuple{F, T}}) where {F, T} = T



is_nothing(::Nothing) = true
is_nothing(::Any) = false

nonscalar_types = unique(filter(!is_nothing, try_get_single_argtype(Base.iterate)))


nonscalar_types = filter(nonscalar_types) do T
    !occursin(".", string(T)) #Quick hack to see if it is in a module that is loaded
end



string_form(T) = string(T)
string_form(T::UnionAll) = string(T.body, " where ", T.var)


function scalarness_function_expr(T)
    str_T = string_form(T)
    parts = split(str_T, "where"; limit=2)
    scalarness_function_expr(parts...)
end


# 1 arg means no where clause
function scalarness_function_expr(type_str::AbstractString)
    type_expr = Meta.parse(type_str)
    
    :(scalarness(::Type{<:$(type_expr)}) = NonScalar())
end

# 2 arg means found a where clause
function scalarness_function_expr(type_str::AbstractString, where_str::AbstractString)
    type_expr = Meta.parse(type_str)
    Meta.parse(string(:(scalarness(::Type{<:$(type_expr)}))) * " where $(where_str) = NonScalar()")
end



for T in nonscalar_types
    try
        eval(scalarness_function_expr(T))
    catch err
        println()
        @show err
        @show T
        @show scalarness_function_expr(T)
    end        
end

@show aslist(Set([1,2]));
dic=Dict(1=>"one",2=>"two")
println(aslist(dic))


