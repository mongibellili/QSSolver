using CodeCosts
f(x::T) where T = convert(T, max(x * 10.0, x / 3))
display(@code_costs f(1.0f0))