
using StaticArrays
using BenchmarkTools


macro testMacro()
  s=0.0
  v1 = @SVector rand(Float64, 100)
  for i = 1:1e+2
    for j = 1:100
        n=rand(1:100)
        s=s+v1[n]
     end
  end
  s
end
function testFunc()
  s=0.0
  v1 = @SVector rand(Float64, 100)
  for i = 1:1e+2
    for j = 1:100
        n=rand(1:100)
        s=s+v1[n]
     end
  end
  s
end

@btime @testMacro()
@btime testFunc()


#display(@simulatorMvector())
#display(simulatorvector())



#=
macro simulatorMvector()
  s=0.0
  v1 = @MVector rand(Float64, 100)
  for i = 1:1e+2
    for j = 1:100
        n=rand(1:100)
        s=s+v1[n]
     end
  end
  s
end
macro simulatorvector()
  s=0.0
  v1 = rand(100)
  for i = 1:1e+2
    for j = 1:100
        n=rand(1:100)
        s=s+v1[n]
     end
  end
  s
end
=#