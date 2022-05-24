using StaticArrays
using BenchmarkTools
#TP7: compare stack-recreate vs heap-update via vectors and svectors and Mvectors
function simulatorSvector()
  #v1 = SVector{:($n),Float64}(undef)
  v1 = @SVector zeros(100)
  for i = 1:1e+3
    for j = 1:100
      v1=setindex(v1,(i - j) * 0.54,j)
    end
  end
  v1
end
function simulatorMvector()
  #v1 = SVector{:($n),Float64}(undef)
  v1 = @MVector zeros(100)
  for i = 1:1e+3
    for j = 1:100
      v1[j] = (i - j) * 0.54
    end
  end
  v1
end
function simulatorvector()
  v1 = Vector{Float64}(undef, 100)
  for i = 1:1e+3
    for j = 1:100
      v1[j] = (i - j) * 0.54
    end
  end
  v1
end
@btime simulatorSvector()
@btime simulatorvector()
@btime simulatorMvector()


#display(simulatorSvector());println()
#display(simulatorvector())

