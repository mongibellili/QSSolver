#= rhs=((schema.args[i].args[2]))
#println(typeof(rhs))
basi = convert(Basic, rhs)#:(2u1+3u2))
#extract jaco components
for j=1:T            # the number of diffequ always coincides with the number of continuous vars?
    coef=diff(basi,u[j])
    push!(jacArr,coef)
end =#
#= using StaticArrays
p = @SVector fill(NaN, 3)


#= f=fill(NaN, 3) =#
display(p);println()
display(typeof(p[1])) =#
#= x=0.0
println(x!==NaN)
#println(x isa NaN)# isa expects a type
y=NaN
println(y!==NaN)
#= println(y==NaN) =#

println(sign(x))
println(typeof(sign(x))) =#
#= using BenchmarkTools
function comp(x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})
    
   # @. x=y+x*3+y*z
     x.=(y.+x.*3).+y.*z

    #println("x= ",x)
end

x = [0.0,2.0,1.0]

y = [1.0,0.0,3.0]

z=[3.5, 4.7, 0.0]

@btime comp(x,y,z)
#@show x =#
i=1
#= s=Symbol("f", i)
println(s) =#
function sd()
    println(i)
end

s=sd
@eval $s()