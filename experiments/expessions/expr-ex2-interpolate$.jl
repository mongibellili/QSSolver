x = :(1 + 2)

a = quote $x end
#display(a)#quote 1 + 2 end
#display(eval(a))#3
b = quote quote $x end end
#display(b)#quote    $(Expr(:quote, quote $(Expr(:$, :x))end))    end
#display(eval(b))#quote 1+2  end
c = quote $(quote $x end) end
#display(c)#quote begin 1 + 2 end   end
#display(eval(c))#3
d = quote quote $$x end end
#display(d)#quote   $(Expr(:quote, quote    $(Expr(:$, :(1 + 2))) end))    end
#display(eval(d))#quote 3 end
