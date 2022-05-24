
using StaticArrays
using BenchmarkTools

mutable struct DataModel
  x::MVector{2,Float64} 
  y::Int

  function DataModel(u::MVector{2,Float64},y::Int)
    DataModel
  end
  function testFun(::Val{1},u::MVector{2,Float64})
    u[1]
  end
end


v1=@MVector [1.0,2.0]
dm = DataModel(v1, 5)
v2=@MVector [1.0,2,0]
display(dm.testFun(Val(1),v1))





