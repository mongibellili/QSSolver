funs = []


 @generated make_fun(::Val{N}, x) where N = funs[N]
#make_fun (generic function with 3 methods)

evil(code) = (push!(funs, code); x->make_fun(Val{length(funs)}(), x))
#evil (generic function with 1 method)

@show evil(:(x*2))(10)