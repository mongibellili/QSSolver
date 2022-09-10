#using ExprTools
using MacroTools
function get_param_name(expr) :: Symbol
    @capture(expr, arg_<:arg_type_) && return arg
    @capture(expr, arg_) && return arg
  end

  function get_args(func_def::Dict)
    arg_dict = Dict{Symbol, Any}()
    arg_list = Vector{Symbol}()
    kwarg_list = Vector{Symbol}()
    for arg in (func_def[:args]...,)
      arg_def = splitarg(arg)
      push!(arg_list, arg_def[1])
      arg_dict[arg_def[1]] = arg_def[3] ? Any : arg_dict[arg_def[1]] = arg_def[2]
    end
    for arg in (func_def[:kwargs]...,)
      arg_def = splitarg(arg)
      push!(kwarg_list, arg_def[1])
      arg_dict[arg_def[1]] = arg_def[3] ? Any : arg_dict[arg_def[1]] = arg_def[2]
    end
    arg_list, kwarg_list, arg_dict
  end


macro log_trace(expr)
    expr.head !== :function && error("Expression is not a function definition!")
    func_def=splitdef(expr)
   # @show func_def
    rtype = :rtype in keys(func_def) ? func_def[:rtype] : Any
   # @show rtype
   # display(def);println()
   args, kwargs, arg_dict = get_args(func_def)
 #  @show args  #[:x, :y, :v]
 #  @show kwargs  #[:bl]
  # @show arg_dict   #Dict{Symbol, Any}(:bl => :Bool, :y => :Float64, :v => :(Vector{T}), :x => :Float64)
  @show func_def[:whereparams]  #(:(T <: Float64),)
  println(func_def[:whereparams]...)  #T <: Float64
  # params = ((param for param in func_def[:whereparams])...,) #(:(T <: Float64),)
   params = ((get_param_name(param) for param in func_def[:whereparams])...,)  #(:T,)
 #  @show params
  #ui8 = BoxedUInt8(zero(UInt8))
   

    #combinedef(def)

end
@log_trace function foo(x::Float64,y::Float64,v::Vector{T};bl::Bool=false )::Float64 where {T<:Float64}
    return y+x
end
