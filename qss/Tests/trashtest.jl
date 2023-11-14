using SymEngine
x=4.0
d=[1.0,1.0,2.0,3.0]
#m=quote 3+5x+d*x end

macro test(m)
Base.remove_linenums!(m)
str=string(m)
@show str
basi = convert(Basic, m)
coef = diff(basi, :x)
@show coef
end

@test 3+5x+d[g]*x