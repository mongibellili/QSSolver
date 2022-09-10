
 macro assert(ex)
    println(ex)
    return :( $ex ? nothing : throw(AssertionError($(string(ex)))) )
end
function assert(ex)
    return :( $ex ? nothing : throw(AssertionError($(string(ex)))) )
end
#display(@macroexpand @assert rf)
#display( @assert 1 == 12.0)#error throw with clear message
#display( @assert 1 == 1.0)#nothing
#display( assert(1 == 12.0))#:(if false nothing else throw(AssertionError("false")) end)
#display( assert(1 == 1.0))#:(if true nothing else throw(AssertionError("true")) end)
#= the expression in the macro gets constructed and placed in the syntax tree where the macro call occurs, 
    the var ex does care if it exists or what its value is...that is runtime's business=#

display(@assert rf)# at parse time rf is printed as is (does not care), at runtime  UndefVarError: rf not defined.