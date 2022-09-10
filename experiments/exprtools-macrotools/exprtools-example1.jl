using ExprTools
ex = :(
           function Base.f(x::T, y::T) where T
               x + y
           end
       )
def = splitdef(ex)
#eval(ex) # LoadError: UndefVarError: f not defined
#display(def)  #outputs: 
                        #= Dict{Symbol, Any} with 5 entries:
                                :args => Any[:(x::T), :(y::T)]
                                :body => quote x + y    end
                                :name => :(Base.f)
                                :head => :function
                                :whereparams => Any[:T] =#
def[:head]= :(=)
#g_expr = combinedef(def)
#eval(g_expr) #LoadError: UndefVarError: f not defined : if i remove the word Base above the error goes away
def[:name] = :g;
g_expr = combinedef(def)
#display((g_expr));println() # :((g(x::T, y::T) where T) = begin x + y end)
display(eval(g_expr));println() #g (generic function with 1 method)
#display(signature(eval(g_expr)))#LoadError: MethodError: no method matching signature(::typeof(g))
display(signature(first(methods(g)))) # Dict{Symbol, Any} with 3 entries: :name => :g  :args => Expr[:(x::T), :(y::T)]    :whereparams => Any[:T] 