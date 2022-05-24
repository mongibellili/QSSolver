using SymEngine

T=4
uarray=[symbols("u$i") for i in 1:T]
display(uarray);println()
s=Symbol(:u3)
display(s)
basi = convert(Basic, s)
positionArray = indexin(basi, uarray)
display(positionArray);println()