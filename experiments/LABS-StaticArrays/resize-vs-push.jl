using StaticArrays
using BenchmarkTools
states=2
initarr = []

for i = 1:states
  push!(initarr, zeros(12))
end
arr1=zeros(12)
savedVars = SVector{states,Array{Float64}}(tuple(initarr...))
initarr2 = []
for i = 1:states
    push!(initarr2, [])
  end
 
  savedVars2 = SVector{states,Array{Float64}}(tuple(initarr2...))


function test(arr::Vector{Float64})
    count = 1
    while count < 10 #&& count < 4
    count += 1
    for k = 1:2
        arr[count]=1.1
    end
    end#end while
    #display(savedVars)
    arr
end

function test2(arr::SVector{2,Array{Float64}})
    count = 1
    while count < 10 #&& count < 4
      count += 1
    for k = 1:2
         arr[k][count]=1.1
     end
    end#end while
    #display(savedVars)
    arr
end
function test3(arr::SVector{2,Array{Float64}})
    count = 1
    while count < 10 #&& count < 4
      count += 1
    for k = 1:2
         push!(arr[k],1.1)
     end
    end#end while
    #display(savedVars)
    arr
end
@btime test(arr1)
@btime test2(savedVars)
@btime test3(savedVars2)