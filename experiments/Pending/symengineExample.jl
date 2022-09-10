
using SymEngine
@vars a 
@vars b
s= Symbol(:u,1)
f1=diff(a + 2(b+2)^2 + 2a + 3(a+1), b)
display(f1)

u=[symbols("u$i") for i in 1:3]
#@vars Symbol(:u,1)   Symbol(:u,2)
#display(typeof(u[1]))
i=2

basi = convert(Basic, :(2u1+3u2+5))
display(length(basi));println()
#= dif=diff(basi,u[i])
display(dif) =#
for i=1:2
global basi=subs(basi,u[i]=>1)
end
display(basi)
#= ev=SymEngine.evalf(basi,3,3)
display(ev) =#