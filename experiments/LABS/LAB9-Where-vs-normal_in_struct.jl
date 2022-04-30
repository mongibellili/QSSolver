
using StaticArrays
using BenchmarkTools


#TP9: compare stack-fetch vs heap-fetch via vectors and svectors in immutable structs
struct unknownDataM
  x::MVector{N,Float64} where {N}  
  y::Int
end


struct DataM
  x::MVector{5,Float64} 
  y::Int
end
struct DataNormal
  x::Vector{Float64} 
  y::Int
end



function testknownM(x::MVector{5,Float64} )
  p1 = DataM(x, 5)
  for i = 1:1e+1 for j = 1:5
      p1.x[j] =p1.x[j]* (-j)  
  end end
end
function testknownM(x::Vector{Float64} )
  p1 = DataM(x, 5)
  for i = 1:1e+1 for j = 1:5
      p1.x[j] =p1.x[j]* (-j)  
  end end
end
function testNormal(x::Vector{Float64} )
  p1 = DataNormal(x, 5)
  for i = 1:1e+1 for j = 1:5
      p1.x[j] =p1.x[j]* (-j)  
  end end
end
function testUnknownDataM(x::MVector{N,Float64} ) where {N}
  p1 = unknownDataM(x, 5)
  for i = 1:1e+1 for j = 1:5
      p1.x[j] =p1.x[j]* (-j)   
  end end
end

#v1 = Vector{Float64}(undef,n)
n=5
v1=rand(n)
v2 =  @MVector rand(n)
#@btime DataNormal(v1, 5)
#@btime DataM(v1, 5)
#@btime unknownDataM(v2, 5)

#@btime testknownM(v1)
#@btime testknownM(v2)
@btime testNormal(v1)
@btime testUnknownDataM(v2)





#= #v1 =  @MVector[1.1,2.2,3.3,4.4,5.5]
n=5
v1 =  @MVector zeros(n)
#= p1 = Data(v1, 5)
p2 = unknownData(v1, 5) =#
#= testWhereEffect(v1)
display(v1);println() =#
#= testNormal(v1)
display(v1);println() =#
@btime testNormal(v1)
@btime testWhereEffect(v1) =#

#display(p2.x)
#= @btime Data(v1, 5)
@btime unknownData(v1, 5) =#





#=
struct DataSvector
  x::SVector{5,Float64} 
  y::Int
end
struct DataMvector
  x::MVector{5,Float64} 
  y::Int
end
struct DataNormal
  x::Vector{Float64}
  y::Int
end
function testSvector()
  s=0.0
  v1 =  @SVector[1.1,2.2,3.3,4.4,5.5]
  p1 = DataSvector(v1, 5)
  for i = 1:1e+3
    for j = 1:5
      s=s+p1.x[j]    
    end
  end
  s
end
function testMvector()
  s=0.0
  v1 =  @MVector[1.1,2.2,3.3,4.4,5.5]
  p1 = DataMvector(v1, 5)
  for i = 1:1e+3
    for j = 1:5
      s=s+p1.x[j]    
    end
  end
  s
end
function testvector()
  s=0.0
  v2 = [1.1,2.2,3.3,4.4,5.5]
  p2 = DataNormal(v2, 5)
  for i = 1:1e+3 for j = 1:5
        s=s+p2.x[j]
  end end
  s
end


@btime testSvector()
@btime testvector()
@btime testMvector()
#display(testSvector());println()
#display(testvector())

=#