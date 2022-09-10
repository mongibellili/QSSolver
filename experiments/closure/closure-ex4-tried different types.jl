using ForwardDiff
function constructF()
    x=0
    for T in (:Int, :Float64)
        @eval function f(x::$T)# if @eval removed: loaderror local variable T cannot be used in closure declaration
                x += 1# the value x does not get stored as in closure ...@eval      
        end
    end
    return f
end

#does not increment
println(constructF())
println(f(0.5))#1.5
println(f(0.5))#1.5
println(f(0.5))#1.5
