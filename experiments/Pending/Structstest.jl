
using StaticArrays
using BenchmarkTools

#=
@generated function funcsimulatorMvector(n)
     m=:(n)
end
using StaticArrays
using BenchmarkTools

mutable struct Data
  x::MVector{N,Float64} where {N}
  y::Int
  function Data(n::Int,y1::Int)
    p=new()
    p.x = @MVector zeros(n)
    p.y=y1
    p
  end


end
mutable struct DataNormal
  x::Vector{Float64}
  y::Int
  function DataNormal(n::Int,y1::Int)
    p=new()
    p.x = Vector{Float64}(undef, n)
    p.y=y1
    p
  end
end

function simulatorMvector(n)

  p = Data(n, 5)
  for j=1:n
  for i=1:100
   p.x[j] = (i-j)*0.54
  end
  end
end
function simulatorvector(states)

  p = DataNormal(states, 5)
  for j=1:states
    for i=1:100
     p.x[j] = (i-j)*0.54
    end
    end
end

@btime simulatorvector(2)
@btime simulatorMvector(2)

=#

struct Data
  x::MVector{N,Float64} where {N}
  y::Int
end
struct DataNormal
  x::Vector{Float64}
  y::Int
end

function simulatorMvector()
  # m=:($n)
  v1 = MVector{2,Float64}(undef)
  p = Data(v1, 5)
  for j = 1:2
    for i = 1:100000
      p.x[j] = (i - j) * 0.54
    end
  end
 #p.y=2

end

function simulatorvector()
  v1 = Vector{Float64}(undef, 2)
  p = DataNormal(v1, 5)
  for j = 1:2
    for i = 1:100000
      p.x[j] = (i - j) * 0.54
    end
  end

end
function simulatorvector()
  v1 = Vector{Float64}(undef, 2)
  p = DataNormal(v1, 5)
  for j = 1:2
    for i = 1:100000
      p.x[j] = (i - j) * 0.54
    end
  end

end

@btime simulatorvector()
@btime simulatorMvector()

#display(@simulatorMvector(2))
#display(@simulatorvector(2))

