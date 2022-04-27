
using StaticArrays
using BenchmarkTools


#TP9: compare stack-fetch vs heap-fetch via vectors and svectors in immutable structs


struct DataMvector
  x::MVector{5,Float64} 
  y::Int
end


function simulatorMvector(x::MVector{5,Float64} )
    for j = 1:5
      n=rand(1:5)
      x[n] = (-j) * 0.54   
    end
end

v1 =  @MVector[1.1,2.2,3.3,4.4,5.5]
p1 = DataMvector(v1, 5)

simulatorMvector(p1.x)
display(v1);println()
display(p1.x)

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
function simulatorSvector()
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
function simulatorMvector()
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
function simulatorvector()
  s=0.0
  v2 = [1.1,2.2,3.3,4.4,5.5]
  p2 = DataNormal(v2, 5)
  for i = 1:1e+3 for j = 1:5
        s=s+p2.x[j]
  end end
  s
end


@btime simulatorSvector()
@btime simulatorvector()
@btime simulatorMvector()
#display(simulatorSvector());println()
#display(simulatorvector())

=#