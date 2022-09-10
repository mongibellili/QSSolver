a=Expr(:$, :(1+2))
#display(dump(a))# head symbol$  args  expression
#display(a)#:($(Expr(:$, :(1 + 2))))
#@show eval(a)#LoadError: syntax: "$" expression outside quote  ; when head=$ causes errors, it can be fixed below (b or c)
b=Meta.quot(a)  #wrap the complainer a  with special expr of head = quote
#@show typeof(b)#Expr
#display(dump(b))#head symbol quote  args  expression=a  =(head=symbol$ args=expression)
#@show b # :( $(Expr(:quote, :($(Expr(:$, :(1 + 2)))))) )
#@show eval(b) #3
c=QuoteNode(a)  #wrap the complainer a  with special object (similar to expr) called QuoteNode that has value=Expr head= symbol$ args=expression
#@show c  # :( $(QuoteNode( :($(Expr(:$, :(1 + 2))))))  )
#display(dump(c))
#@show typeof(c)#QuoteNode
d=eval(c)# evaluates to a
#d= :($(Expr(:$, :(1 + 2))))# exactly a
#@show typeof(d)#Expr
#@show (d==a)#true