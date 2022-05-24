using ForwardDiff
function constructF(g::Function, n::Int)

    v = [g]
    for i in 2:n
        gi(x) = ForwardDiff.derivative(v[i - 1], x)
        println(gi(1.2))
        v = vcat(v, gi)
    end

    function f(x::T) where{T<:Number}
        acc = zero(T)
        for i in 1:n
            acc += v[i](x)^i
        end
        acc
    end

    return f
end

g(x)=x*x
h=constructF(g,2)
display(h(2.1))