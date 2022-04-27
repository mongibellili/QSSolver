using StaticArrays
using BenchmarkTools
#TP10: macro-compare stack-recreate vs heap-update via vectors and svectors and Mvectors and random
macro simulatorSvector()
  #v1 = SVector{:($n),Float64}(undef)
  v1 = @SVector rand(Float64, 100)
  for i = 1:1e+3
    for j = 1:100
      n=rand(1:100)
      v1=setindex(v1,(i - j) * 0.54,n)
    end
  end
  v1
end
macro simulatorMvector()
  #v1 = SVector{:($n),Float64}(undef)
  v1 = @MVector rand(Float64, 100)
  for i = 1:1e+3
    for j = 1:100
      n=rand(1:100)
      v1[n] = (i - j) * 0.54
    end
  end
  v1
end
macro simulatorvector()
  v1 = rand(100)
  for i = 1:1e+3
    for j = 1:100
      n=rand(1:100)
      v1[n] = (i - j) * 0.54
    end
  end
  v1
end
@btime @simulatorSvector()
@btime @simulatorvector()
@btime @simulatorMvector()


#display(@simulatorSvector());println()
#display(@simulatorvector());println()
display(@simulatorMvector())

