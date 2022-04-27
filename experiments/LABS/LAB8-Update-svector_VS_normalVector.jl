using StaticArrays
using BenchmarkTools
#TP9: vectors created out
function simulatorSvector(v1::SVector{10,Float64} ) 
  for i = 1:1e+3
    for j = 1:10
      n=rand(1:10)
      v1=setindex(v1,(i - j) * 0.54,n)
    end
  end
  v1
end
function simulatorMvector(v1::MVector{100,Float64} ) 
  for i = 1:1e+3
    for j = 1:10
      n=rand(1:10)
      v1[n] = (i - j) * 0.54
    end
  end
  v1
end
function simulatorvector(v1::Vector{Float64})
  for i = 1:1e+3
    for j = 1:10
      n=rand(1:10)
      v1[n] = (i - j) * 0.54
    end
  end
  v1
end
v1 = @SVector rand(Float64, 10)
v2 = @MVector rand(Float64, 10)
v3 = rand(10)
@btime simulatorSvector(v1)
@btime simulatorvector(v3)
@btime simulatorMvector(v2)


#display(simulatorSvector());println()
#display(simulatorvector())

