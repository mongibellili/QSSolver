foo = Dict(:a=>1, :b=>2, :c=>3)
println("**********tuto1*********")
function m(; kwargs...)   # if delete 3 dots then LoadError: UndefKeywordError: keyword argument kwargs not assigned
    display(kwargs);println()  #pairs(::NamedTuple) with  entries: :symbol => value etc...
   #=  println(kwargs.a)
    println(a) =#
    display((;kwargs...))  # change to named tuple   #(keyname = value)
end
#display(m(foo))#LoadError: MethodError: no method matching m(::Dict{Symbol, Int64})
#display(m(foo...))#LoadError: MethodError: no method matching m(::Pair{Symbol, Int64}, ::Pair{Symbol, Int64}, ::Pair{Symbol, Int64})
#m(;foo)     #outputs:
                    #= pairs(::NamedTuple) with 1 entry:
                    :foo => Dict(:a=>1, :b=>2, :c=>3)
                    (foo = Dict(:a => 1, :b => 2, :c => 3),) =#
# m(;foo...)     #outputs:
                    #=  pairs(::NamedTuple) with 3 entries:
                    :a => 1 :b => 2  :c => 3
                    (a = 1, b = 2, c = 3) =#

println("**********tuto2*********")
function m2(; kwargs...)
    println(kwargs)
    println(kwargs[:a])
end
#display(m2(;foo))#LoadError: type NamedTuple has no field a
#display(m2(;foo...))  Base.Pairs(:a => 1, :b => 2, :c => 3)       1
println("**********tuto3*********")
m3(;kwargs...) = dump(kwargs)
#m3(;foo...)  #outputs:
                #= Base.Pairs{Symbol, Int64, Tuple{Symbol, Symbol, Symbol}, NamedTuple{(:a, :b, :c), Tuple{Int64, Int64, Int64}}}
                data: NamedTuple{(:a, :b, :c), Tuple{Int64, Int64, Int64}}
                a: Int64 1  b: Int64 2   c: Int64 3
                itr: Tuple{Symbol, Symbol, Symbol}
                1: Symbol a 2: Symbol b 3: Symbol c =#
m4(;kwargs...) = kwargs.data   #Warning: use values(kwargs) and keys(kwargs) instead of kwargs.data and kwargs.itr
m5(;kwargs...) = values(kwargs)
m6(;kwargs...) = keys(kwargs)
#display(m5(;foo...)) #(a = 1, b = 2, c = 3)
#display(m6(;foo...)) #(:a, :b, :c)