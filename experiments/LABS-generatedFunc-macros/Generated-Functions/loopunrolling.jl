using BenchmarkTools
function avg2(vals::Vector{T})where {T}
    sum = vals[1]
    for i in 2:length(vals)
          sum += vals[i]
    end
    sum/length(vals)
end

struct Vector2{N,T}
    vals::Vector{T}
end
@generated function avgg(els::Vector2{N,T}) where {N,T}
    code = :(els.vals[1])
    for i=2:N
        code = :($code + els.vals[$i])
    end
    :(($code)/$N)
end
s = Vector2{10,Int64}([1,2,3,4,9,12,2,5,3,7])

@btime avg2($s.vals)
@btime avgg($s)