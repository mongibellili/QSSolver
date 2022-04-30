
using StaticArrays
########################################Normal vectors##########################################
#= l=Array[] #Vector{Array}
m=Float64[]# Vector{Float64}
n=Array{Float64}(undef,0)# Vector{Float64}
n1=rand(5)
o= Array{Float64}(undef,1)# Vector{Float64}
p= Array{Float64,2}(undef,1,2)# Matrix{Float64}

#r=Array{Float64}[]#Vector{Array{Float64}} but it has 0 elements r[1] is illegal
#s=[1]
#t=zeros(Float64, 2,2)

#=r=Vector{Array{Float64}}(undef, 1) # it has 1 Array element
display(r);println()
r[1]=Array{Float64}[]
push!(r[1],0)
@show typeof(r[1])=#

#@show typeof(l);println()
#@show typeof(m);println()
#@show typeof(n1);println()
#@show typeof(o);println()
#@show typeof(p);println()
#=
@show typeof(r);println()
@show typeof(s);println()
@show typeof(t);println()
#display(s)
=# =#


########################################Svectors##########################################
#= 
v0=SA[1,2,3]
#display(v0[1]);println()
v1 = SVector(1,2,3)
#display(v1);println()
v2 = @SVector [1,2,3]
#display(typeof(v1));println()
#display(v2[1]);println()
v3=SVector{3,Int}(1,2,3)
#display(v3[1]);println()
vect=(1,2,3)
v4=SVector{3,Int}(vect)
#display(v4[1]);println()
v5 = @SVector zeros(3)
#display(v5[1]);println()
v5=setindex(v5,10,2)
#display(v5);println()
v6 = @SVector randn(Float64, 40)
#display(size(v1));println()
struct Data
    x::SVector{N,Float64} where {N}
    y::Int
end
n=2
v1 = @SVector zeros(n)
p = Data(v1, 5)
#display(p.x);println() =#
#= states=2
arr=[]
#for i=1:states

savedVars=SVector{2,Array{Float64}}([],[])
#display(typeof(savedVars))
push!(savedVars[1],0.3)
display(savedVars) =#
#= n=2
v2 =  @MVector randn(Float64,n)
display(v2);println() =#

#SVector{nColumns,Array{Int}}
v1 = tuple( [1,2,3]...)
v2 = tuple(  [1,3]...)
v3 =  tuple( [1,3,4,5]...)
arr=[]
push!(arr,v1)
push!(arr,v2)
push!(arr,v3)
#arr=[v1,v2,v22]
#display(arr)
#push!(arr,v1)
#= for i = 1:states 
    push!(arr,[])        
end =#
v4= SVector{3, SArray{Tuple, Int64, 1}}
v5=(tuple(arr...)) 
#v3 = @SVector [v1,v2,v22]
#display(typeof(v5));println()
display(v5[1][2]);println()
#dep=SVector{3, SArray{S, Int64, 1} where S<:Tuple}
#dep=setindex(dep,[1,2],1) error
#display(dep[1]);println()  error
######################################## MVector #########################################
#= v1 = MVector(1,2,3)
#display(v1);println()
display(typeof(v1));println()
v2 = @MVector [1,2,3]
#display(v2);println()
v3=MVector{3,Int}(1,2,3)
#display(v3);println()
vect=(1,2,3) #or [1,2,3]
v4=MVector{3,Int}(vect)
#display(v4);println()
v5 = @MVector zeros(3) #Float64
#display(v5);println()
v6 = @MVector rand(Float64, 40)
#display(typeof(v6));println()
#display(v6) =#
#= v1 = MVector{2,Float64}(undef)
display(typeof(v1));println()
v1[1]=2.0
display(v1);println() =#

########################################SMatrix#########################################
#2×2 SMatrix{2, 2, Int64, 4} 
#=  m1 = SMatrix{2,2}(1, 2, 3, 4)
display(m1);println() 
m2 = @SMatrix[1 2;3 4]
display(m2);println() =#
#  error: no constructor m1 = SMatrix(1,2,3,4);m1 = @SMatrix(1,2,3,4)
#=  m3 = SMatrix{2,2}([1 3 ; 2 4])
display(m3);println() 
m3 = transpose(m3)
display(m3);println()  =#