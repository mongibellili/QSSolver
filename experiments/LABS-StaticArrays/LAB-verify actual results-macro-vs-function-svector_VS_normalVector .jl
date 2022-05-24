using StaticArrays
using BenchmarkTools
#verify same results : function vs macro results -compare stack-recreate vs heap-update via vectors and svectors and Mvectors and random
macro simulatorSvector()
  #v1 = SVector{:($n),Float64}(undef)
  v1 =  @SVector[1.1,2.2,3.3,4.4,5.5]
  for i = 1:1e+1
    for j = 1:5 
      v1=setindex(v1,(i - j) * 0.54,j)
    end
  end
  v1
end
macro simulatorMvector()
  #v1 = SVector{:($n),Float64}(undef)
  v1 =  @MVector[1.1,2.2,3.3,4.4,5.5]
  for i = 1:1e+1
    for j = 1:5
        v1[j] = (i - j) * 0.54
    end
  end
  v1
end
macro simulatorvector()
  v1 = [1.1,2.2,3.3,4.4,5.5]
  for i = 1:1e+1
    for j = 1:5
      v1[j] = (i - j) * 0.54
    end
  end
  v1
end
function simulatorSvector()
  #v1 = SVector{:($n),Float64}(undef)
  v1 =  @SVector[1.1,2.2,3.3,4.4,5.5]
  for i = 1:1e+1
    for j = 1:5 
      v1=setindex(v1,(i - j) * 0.54,j)
    end
  end
  v1
end
function simulatorMvector()
  #v1 = SVector{:($n),Float64}(undef)
  v1 =  @MVector[1.1,2.2,3.3,4.4,5.5]
  for i = 1:1e+1
    for j = 1:5
        v1[j] = (i - j) * 0.54
    end
  end
  v1
end
function simulatorvector()
  v1 = [1.1,2.2,3.3,4.4,5.5]
  for i = 1:1e+1
    for j = 1:5
      v1[j] = (i - j) * 0.54
    end
  end
  v1
end
#@btime @simulatorSvector()
#@btime @simulatorvector()
#@btime @simulatorMvector()


display(simulatorSvector());println()
display(simulatorvector());println()
display(simulatorMvector())

display(@simulatorSvector());println()
display(@simulatorvector());println()
display(@simulatorMvector())