macro parse_run(arg)
    name=54
    println("at Parse time, The argument is: ", name)
    n=5
    #n=n+arg
    println("new n= ", n)
    :(println("at Runtime?, The argument is: ", ($name)))
     return esc(:(println("at Runtime, The argument is: ", ($arg))))
    
end

#= function parse_run(arg)
    println("from a function at Parse time, The argument is: ", arg)
    return :(println("from a function at Runtime, The argument is: ", $(esc(arg))))
end =#

function test_macro(name::String)
    
   res= @parse_run name
   # parse_run(name)
   println("res=  ",res)
end

name="mongi"
#display(@macroexpand @parse_run name)

test_macro("Amjad")