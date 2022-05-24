using StaticArrays
using BenchmarkTools


  #v1 = MVector{:($n),Float64}(undef)
  #v1 =  @MVector zeros(n)
   # v1=@MVector rand(n)
function testUnknownSize(n)
  v1 = MVector{n,Float64}(undef)
  for j = 1:n
    for i = 1:10
      v1[j] = (i - j) * 0.54
    end
  end
  v1
end
function testknownSize(n)
  #v1 = MVector{2,Float64}(undef)
  v1=@MVector rand(2)
  for j = 1:2
    for i = 1:10
      v1[j] = (i - j) * 0.54
    end
  end
  v1
end
function testNormalVector(n)
  #v1 = Vector{Float64}(undef,2)
  v1 = rand(2)
  for j = 1:2
    for i = 1:10
      v1[j] = (i - j) * 0.54
    end
  end
  v1
end


@btime testUnknownSize(2)
@btime testknownSize(2)
@btime testNormalVector(2)


#= display(@testUnknownSize(2));println()
display(@testUnknownSize(2));println()
display(testknownSize(2)) =#
#display(@macroexpand @testUnknownSize(2))



#= macro testUnknownSize(n)
  #v1 = MVector{:($n),Float64}(undef)
  v1 = MVector{n,Float64}(undef)
 # p = Data(v1, 5)
  for j = 1:n
    for i = 1:10
      #p.x[j] = (i - j) * 0.54
      v1[j] = (i - j) * 0.54
    end
  end
  v1
end =#