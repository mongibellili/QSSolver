using StaticArrays
using BenchmarkTools




struct Data
  x::MVector{N,Float64} where {N}
  y::Int
end
macro simulatorMvector(n)
    #v1 = MVector{funcsimulatorMvector(n),Float64}(undef)#did not improve performance
    v1 = MVector{:($n),Float64}(undef)
    p = Data(v1, 5)
    for j = 1:n
      for i = 1:100
        p.x[j] = (i - j) * 0.54
      end
    end
    #p.x
end
function simulatorMvector(n)
  #v1 = MVector{funcsimulatorMvector(n),Float64}(undef)#did not improve performance
  v1 = MVector{n,Float64}(undef)
  p = Data(v1, 5)
  for j = 1:n
    for i = 1:100
      p.x[j] = (i - j) * 0.54
    end
  end
  p.x
end
#@btime @simulatorMvector(2)
#@btime simulatorMvector(2)
#display(@simulatorMvector(2));println()
#display(simulatorMvector(2))
display(@macroexpand @simulatorMvector(2))