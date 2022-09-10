
using InteractiveUtils

replace_sin(x::Symbol) = x == :sin ? :cos : x
replace_sin(e::Expr) = Expr(e.head, map(replace_sin, e.args)...)
replace_sin(u) = u

#= macro replace_sin(ex)
	replace_sin(esc(ex))
end =#

function outerreplace_sin(ex)
	replace_sin((ex))
end

#display(@replace_sin(cosp1(x) = 1 + sin(x)));println()
#display(@macroexpand @replace_sin(cosp1(x) = 1 + sin(x)))
#display(@code_lowered @replace_sin(cosp1(x) = 1 + sin(x)))
myexpr=:(cosp1(x) = 1 + sin(x))
res=outerreplace_sin(myexpr)
display(eval(res))
display(cosp1(0.009))

#= function cosp2(x)
	@replace_sin 2 + sin(x)
end
display(cosp2(1)) =#