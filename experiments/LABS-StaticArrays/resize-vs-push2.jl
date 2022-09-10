using StaticArrays
using BenchmarkTools
function test1()
  states=2
  arr = Vector{Array{Float64}}(undef, 2)
  for i = 1:states
      arr[i]=Float64[]
    end
    count = 0
    while count < 2#50
      count += 1
    for k = 1:2
         push!(arr[k],1.1)
     end
    end
    arr
end#known size normal vect push
function test2()
  states=2
  initarr2 = Vector{Array{Float64}}()
  for i = 1:states
    push!(initarr2, Float64[])
  end
  arr = SVector{states,Array{Float64}}(tuple(initarr2...))
    count = 0
    while count < 2#50
      count += 1
    for k = 1:2
         push!(arr[k],1.1)
     end
    end
    arr
end#svector
function test3()
  states=2
  arr = Vector{Array{Float64}}()
  for i = 1:states
    push!(arr, Float64[])
  end
    count = 0
    while count < 2#50 
      count += 1
    for k = 1:2
         push!(arr[k],1.1)
     end
    end
    arr
end#normal vect push
function test4()
  states=2
  arr = Vector{MVector{30,Float64}}(undef, 2)
   for i = 1:states
      arr[i]=zeros(2)
    end
   count = 0
    while count < 2#50 
      count += 1
    for k = 1:2
         arr[k][count]=1.1
     end
    end
    arr
end#unpractical
function test5()
  states=2
#=   arr = Vector{Vector{Float64}}(undef, 2)
   for i = 1:states
      arr[i]=zeros(15)
    end =#
  arr = Vector{Array{Float64}}()
  for i = 1:states
    push!(arr, zeros(2))
  end
   count = 0
   len=length(arr[1])
    while count < 2#50 
      count += 1
    for k = 1:2
      if len<count
        for i=1:2
          resize!(arr[i],count*2)
        end
        len=count*2
      end
         arr[k][count]=1.1
     end
    end
    arr
end
function test6()
  states=2
  arr = Vector{Vector{Float64}}(undef, 2)
   for i = 1:states
      arr[i]=zeros(150)
    end
   count = 0
   len=length(arr[1])
    while count < 250#50 
      count += 1
    for k = 1:2
      if len<count
        for i=1:2
          resize!(arr[i],count*2)
        end
        len=count*2
      end
         arr[k][count]=1.1
     end
    end
    arr
end
function test7()
  states=2
  arr = Vector{Vector{Float64}}(undef, 2)
   for i = 1:states
      arr[i]=Vector{Float64}(undef, 250)
    end
   count = 0
   len=length(arr[1])
    while count < 250#50 
      count += 1
    for k = 1:2
      if len<count
        for i=1:2
          resize!(arr[i],count*2)
        end
        len=count*2
      end
         arr[k][count]=1.1
     end
    end
    arr
end


#= o= Array{Float64}(undef,5)
display(o);println()

r=Vector{Array{Float64}}(undef, 2) # it has 1 Array element
#display(r);println()#2-element Vector{Array{Float64}}: undef undef
push!(r,[1.0]) =#

#= display(test1());println()
display(test2());println()
display(test3());println()
display(test5());println()
display(test6());println() =#
#= @show test1()
@show test2()
@show test3()
@show test5()
@show test6() =#



#= @btime test1()
@btime test2()
@btime test3()
@btime test5()=#
#@btime test6() 
@btime test7()