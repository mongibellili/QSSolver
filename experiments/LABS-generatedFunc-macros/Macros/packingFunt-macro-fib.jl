

function memoit(f::Function,p)
    if !isdefined(Main,:my_memoit_cache)
        global my_memoit_cache =Dict{Function,Dict{Any,Any}}()
    end
    cache = haskey(my_memoit_cache,f) ? my_memoit_cache[f] : my_memoit_cache[f]=Dict()
    haskey(cache,p) ? cache[p] : cache[p] = f(p)
end

macro memo(e)
    (!(typeof(e) <: Expr) || !(e.head == :call)) &&
    error("Wrong @memo params - required a function call")
    dump(e)
#=     return quote
        memoit($(e.args[1]),$(esc(e.args[2])))
    end =#
    (quote
        memoit($(e.args[1]),$(esc(e.args[2])))
    end)
end

#= function fib3(n)
    n <= 2 ? 1 : (@memo fib3(n-1)) + (@memo fib3(n-2))
end =#
function test(str::Int)
    return @memo test2(str)
end
function test2(str::Int)
    return str*(str-1)
end
#display(fib3(6))
display( test(3))