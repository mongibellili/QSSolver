using StaticArrays
#v1 = @SVector randn(Float64, 40)
#display(v1)
#A = rand(3);
#display(A)
#= using Base.Cartesian
display(@macroexpand @nloops 2 i A begin
    r += @nref 2 A i
end) =#
#= a=Float64(0)

if a==0
    println("a is zero")
end =#

#= struct DataSvector
    x::SVector{N,Float64} where {N} 
    y::Int
  end

n=5
v1 = @SVector zeros(n) 
p = DataSvector(v1, 5)  
display(p) =#
#= arr=[[],[],[]]
#= arr=[]
push!(arr,[]) =#
t=tuple(arr...)
display(t) =#
display(sign(5))