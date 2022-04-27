
using StaticArrays
using BenchmarkTools



#TP5: taking randomness outside
function simulatorSvector(v1::SVector{100,Float64} ) 
  s=0.0
  for i = 1:1e+3
    for j = 1:100
      n=rand(1:100)
      s=s+v1[n]
     end
  end
  s
end
function simulatorMvector(v1::MVector{100,Float64} ) 
  s=0.0
  for i = 1:1e+3
    for j = 1:100
      n=rand(1:100)
      s=s+v1[n]
     end
  end
  s
end
function simulatorvector(v1::Vector{Float64})
  s=0.0
  for i = 1:1e+3
    for j = 1:100
      n=rand(1:100)
      s=s+v1[n]
     end
  end
  s
end
v1 = @SVector rand(Float64, 100)
v2 = @MVector rand(Float64, 100)
v3 = rand(100)
@btime simulatorSvector(v1)
@btime simulatorvector(v3)
@btime simulatorMvector(v2)

#display(simulatorMvector(2))
#display(simulatorvector(2))

