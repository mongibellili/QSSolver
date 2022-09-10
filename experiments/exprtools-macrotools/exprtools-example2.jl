using ExprTools

 
    expr = :(x->4*x)
    def=splitdef(expr)
  #  display(def);println()
  #  push!(def,:name => :g)
    def[:name] = :g
    def[:head]= :(=)
    name=Meta.quot(get(def,:name,Symbol("<anon>")))
    def[:body]=quote
        println("entering ", $name, $(args_tuple_expr(def)))
        $(def[:body])
    end
   
   # def[:name] = :g
  
   #display(def);println()
   a= combinedef(def)
   display(a)


