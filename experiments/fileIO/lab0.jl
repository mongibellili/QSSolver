
open("diffequFunc.jl", "w") do io
#open("./experiments/fileIO/diffequFunc.jl", "w") do io
    println(io, "function f1(q::Vector{Float64})")
    println(io, "q[1]+3.0")
    println(io, "end")

end