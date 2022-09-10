#aslist(x)=[x]
aslist(x::AbstractArray)=x
aslist(x::AbstractString)=x
#= println(aslist([1 2 3]))
println(aslist([1,2,3]))
println(aslist(1))
println(aslist("1")) =#
try_get_single_argtype(::Type{Tuple{F, T}}) where {F, T} = T
#for mm in 
   # mm=methods(Base.iterate).ms[1]#[1:5] #methods(Base.iterate)

#println(methods(Base.iterate))
#mm = first(methods(Base.iterate))

mm=methods(Base.iterate)[29]
    sig=mm.sig
    @show sig
    rettype=try_get_single_argtype(sig)
    @show rettype
    @show eltype(rettype)
    @show typeof(rettype)
    aslist(x::rettype)=x

    aslist(x)=[x]
    println(aslist(rettype))
    println(aslist("rettype"))
#end



