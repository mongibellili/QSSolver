using InteractiveUtils
function funny3(n::Int)
    println(" you called f ",n)
end

index=3
s=Symbol("funny",3)
#= display(s);println()
display(@eval $s(4)) =#
#= f = getfield(Main, s)
display(f)
f(10) =#

#---------------------------------------------
for i = 1:3
    f = Symbol(:h_,i)
    #display(@eval $f);println()
    # display(@macroexpand @eval $f(x) = $i*x)
    @eval $f(x) = $i*x
end
f = getfield(Main,Symbol(:h_,2))
display(getfield(Main,Symbol(:h_,2))(3.0));println()
