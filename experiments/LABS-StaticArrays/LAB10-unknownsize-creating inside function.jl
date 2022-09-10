using StaticArrays
using BenchmarkTools
using InteractiveUtils

  #v1 = MVector{:($n),Float64}(undef)
  #v1 =  @MVector zeros(n)
   # v1=@MVector rand(n)
function testUnknownSize(n::Int,q::Vector{Float64})
  df=3.23
  #v1 = MVector{n,Float64}(undef)
  v1=@MVector zeros(n)
  #display(typeof(v1));println()
    for j = 1:n
      v1[j] =abs(q[j]/df)
    end
  
  v1
end
#= function testknownSize(n)
  #v1 = MVector{2,Float64}(undef)
  v1=@MVector zeros(n)
  for j = 1:n
    
      v1[j] = ( j) * 0.54
    
  end
  v1
end =#
function testNormalVector(n::Int,q::Vector{Float64})
  #v1 = Vector{Float64}(undef,2)
  df=3.23
  v1 = zeros(n)
  for j = 1:n
    v1[j] =abs(q[j]/df)
       
    
  end
  v1
end

quantum=[1.1,2.1,3.1,1.0,0.4]
@btime testUnknownSize(5,quantum)
#@btime testknownSize(5,quantum)
@btime testNormalVector(5,quantum)



#display(@code_typed testUnknownSize(2,quantum))
#display(@code_typed testknownSize(2))
#display(@code_typed testNormalVector(2,quantum))


#display(@code_warntype testUnknownSize(2,quantum))
#display(@code_warntype testNormalVector(2,quantum))


#display( testUnknownSize(2,quantum))


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