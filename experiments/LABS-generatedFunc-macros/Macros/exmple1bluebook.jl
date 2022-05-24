#= macro containervar(container,elemt)
    esc(quote
        ($Symbol(container,elemt))=$container[$elemt]
    end)
end =#

macro containervar(container,elemt)
    return (:($(Symbol(container,elemt))=$container[$elemt]))
    
end

display(@macroexpand @containervar lette 2)
