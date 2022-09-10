using MacroTools: prewalk, postwalk, @capture

#= macro foo(ex)
  postwalk(ex) do x
    @capture(x, some_pattern) || return x
    return new_x
  end
end =#

macro foo2(ex)
  postwalk(ex) do x
    @capture(x,ff_(args__)) || return x
    return :($ff(5,$(args...)) ) 
  end
end

ex = quote
  x = f(y, g(z))
  return h(x)
end
# we want to insert an extra argument 5 into all functions 
#= res= postwalk(x -> @capture(x,ff_(args__)) ? :($ff(5,$(args...)) ) : x, ex)
@show res  =#
res=@foo2 quote 
  x = f(y, g(z))
  return h(x)
end
@show res