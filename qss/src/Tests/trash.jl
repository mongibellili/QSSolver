
using BenchmarkTools
function test(t::Float64)
iszero(t)
end
function test2(t::Float64)
  t==0.0
end
t=2.33

@btime test2(t)
@btime test(t)
