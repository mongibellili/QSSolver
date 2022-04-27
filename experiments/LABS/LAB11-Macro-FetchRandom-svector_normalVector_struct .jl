
using StaticArrays
using BenchmarkTools



#TP6: macro-compare stack-fetch vs heap-fetch via vectors  svectors Mvectors and random and random access
macro simulatorSvector()
  # m=:($n)
  #v1 = SVector{:($n),Float64}(undef)
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
macro simulatorMvector()
  # m=:($n)
  #v1 = SVector{:($n),Float64}(undef)
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
  #v1 = Vector{Float64}(undef, 10)
  v1 = rand(100)
  for i = 1:1e+2
    for j = 1:100
        n=rand(1:100)
        s=s+v1[n]
     end
  end
  s
end
@btime @simulatorSvector()
@btime @simulatorvector()
@btime @simulatorMvector()

display(@simulatorMvector())
#display(simulatorvector())

