using ForwardDiff
function constructF(g::Function, n::Int)

    v = [g]
    for i in 2:n
        gi(x) = ForwardDiff.derivative(v[i - 1], x)
        v = vcat(v, gi)
    end
#= println(v[1](2))
println(v[2](2))
println(v[3](2)) =#
    for T in (:Int, :Float64)
        @eval function f(x::$T)# if @eval removed: loaderror local variable T cannot be used in closure declaration
            acc = $T(0)
            for i in 1:$n
                acc += $v[i](x)^i
                @show acc
            end
            acc
        end
    end

    return f
end

g(x)=x*x
println(constructF(g,3))
println(f(0.5))
