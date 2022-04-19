
#=l=Array[] #Vector{Array}
m=Float64[]# Vector{Float64}
n=Array{Float64}(undef,0)# Vector{Float64}
o= Array{Float64}(undef,1)# Vector{Float64}
p= Array{Float64,2}(undef,1,2)# Matrix{Float64}
=#
#r=Array{Float64}[]#Vector{Array{Float64}} but it has 0 elements r[1] is illegal
#s=[1]
#t=zeros(Float64, 2,2)

r=Vector{Array{Float64}}(undef, 1) # it has 1 Array element
display(r);println()
r[1]=Array{Float64}[]
push!(r[1],0)
@show typeof(r[1])

#=@show typeof(l);println()
@show typeof(m);println()
@show typeof(n);println()
@show typeof(o);println()
@show typeof(p);println()
@show typeof(r);println()
@show typeof(s);println()
@show typeof(t);println()
#display(s)
=#
