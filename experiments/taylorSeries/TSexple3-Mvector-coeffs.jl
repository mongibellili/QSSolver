using TaylorSeries
using StaticArrays
#use_show_default(true)
#q=Taylor1(@MVector[2.0,1.3,1.0],4)
q=Taylor1([2.0,1.3,1.0],4)
println(isbits(q))
println((q.coeffs))