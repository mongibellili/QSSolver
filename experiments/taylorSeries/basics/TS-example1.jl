using TaylorSeries
#displayBigO(false)
#use_show_default(true)
l=1
h=4
println("----------------Tutorial -1: diff cosine--------------")
#= t=Taylor1(2)
a=Taylor1([1.0,2.0,5.0])
b= Taylor1([2.2,2.8,3.0])
x=[a,b]

aa=b(t+0.5)
@show aa
a1=evaluate(aa,t+0.5)
@show a1
aaa=aa(t+0.5)
@show aaa
 =#



#= display(t+0.5);println()
@show b
@show (b(t+0.5))#4.35 + 5.8 t + 3.0 t² + 𝒪(t³)
@show b(0.5)
c=differentiate(b)
@show c
@show c(0.5)
d=differentiate(c)
@show d
@show d(0.5) =#
#= function test(x::Vector{Taylor1{Float64}},t::Taylor1{Float64})
    println(x)
    x[2]+t
end

@show test(x,t)
@show test(x,t+0.5)
@show test(x(t+0.5),t+0.5)

@show x[1] =#





#= c=exp(b)
@show c
@show differentiate(c) =#
#println("----------------Tutorial 0: low order different coeffs--------------")
#= n = Taylor1(l)          #@show n            #n =  1.0 t
p = Taylor1([0.0],l)      #@show p            #p =  0
r = Taylor1([1],l)      #@show r            #r =  1
s = Taylor1([0,1],l)    #@show s            #s =  1 t
u = Taylor1([1,0],l)    #@show u            #u =  1
v = Taylor1([1,1],l)    #@show v            #v =  1 + 1 t
w = Taylor1([1,0,0],l)  #@show w            #w =  1 
x = Taylor1([1,0,1],l)  #@show x            #x =  1
y = Taylor1([1,1,0],l)  #@show y            #y =  1 + 1 t
z = Taylor1([1,1,1],l)  #@show z            #z =  1 + 1 t =#
zz = Taylor1([0,0,1],l)  #@show z            #z =  1 + 1 t
#@show exp(zz)
println("----------------Tutorial 1: high order different coeffs--------------")
#= nh = Taylor1(h)          #@show nh         #nh =  1.0 t
ph = Taylor1([0],h)      #@show ph         #ph =  0
rh = Taylor1([1],h)      #@show rh         #rh =  1
sh = Taylor1([0,1],h)    #@show sh         #sh =  1 t
uh = Taylor1([1,0],h)    #@show uh         #uh =  1 
vh = Taylor1([1,1],h)    #@show vh         #vh =  1 + 1 t
wh = Taylor1([1,0,0],h)  #@show wh         #wh =  1
xh = Taylor1([1,0,1],h)  #@show xh         #xh=   1 + 1 t²
yh = Taylor1([1,1,0],h)  #@show yh         #yh =  1 + 1 t
zh = Taylor1([1,1,1],h)  #@show zh         #zh =  1 + 1 t + 1 t² =#
#= zzh = Taylor1([0,0,1],h)  #@show zh         #zh =  1 + 1 t + 1 t²
@show (zzh)
@show sqrt(zzh)
zh = Taylor1([1,1,1],h)  
vh = Taylor1([1,1],h) 
@show sqrt(vh) =#
#=println("----------------Tutorial 2: low order Usage --------------")

@show exp(n)           #n =1.0 t          exp(n) =  1.0 + 1.0 t
@show exp(p)           #p =  0            exp(p) =  1.0
@show exp(r)           #r =  1            exp(r) =  2.718281828459045
@show exp(s)           #s =  1 t          exp(s) =  1.0 + 1.0 t
@show exp(u)           #u =  1            exp(u) =  2.718281828459045
@show exp(v)           #v =  1 + 1 t      exp(v) =  2.718281828459045 + 2.718281828459045 t 
@show exp(w)           #w =  1            exp(w) =  2.718281828459045
@show exp(x)           #x =  1            exp(x) =  2.718281828459045
@show exp(y)           #y =  1 + 1 t      exp(y) =  2.718281828459045 + 2.718281828459045 t    
@show exp(z)           #z =  1 + 1 t      exp(z) =  2.718281828459045 + 2.718281828459045 t =#
#=println("----------------Tutorial 3: high order Usage--------------")

@show exp(nh)        #nh= 1.0 t           exp(nh) =  1.0 + 1.0 t + 0.5 t² + 0.16666666666666666 t³ + 0.041666666666666664 t⁴
@show exp(ph)        #ph= 0               exp(ph) =  1.0
@show exp(rh)        #rh= 1               exp(rh) =  2.718281828459045 
@show exp(sh)        #sh= 1 t             exp(sh) =  1.0 + 1.0 t + 0.5 t² + 0.16666666666666666 t³ + 0.041666666666666664 t⁴
@show exp(uh)        #uh= 1               exp(uh) =  2.718281828459045 
@show exp(vh)        #vh= 1 + 1 t         exp(vh) =  2.718281828459045 + 2.718281828459045 t + 1.3591409142295225 t² + 0.45304697140984085 t³ + 0.11326174285246021 t⁴
@show exp(wh)        #wh= 1               exp(wh) =  2.718281828459045
@show exp(xh)        #xh= 1 + 1 t²        exp(xh) =  2.718281828459045 + 2.718281828459045 t² + 1.3591409142295225 t⁴
@show exp(yh)        #yh= 1 + 1 t         exp(yh) =  2.718281828459045 + 2.718281828459045 t + 1.3591409142295225 t² + 0.45304697140984085 t³ + 0.11326174285246021 t⁴         
@show exp(zh)        #zh= 1 + 1 t + 1 t²  exp(zh) =  2.718281828459045 + 2.718281828459045 t + 4.077422742688568 t² + 3.171328799868886 t³ + 2.8315435713115056 t⁴
=#
#= println("----------------Tutorial 4: low order Usage + 1 --------------")

@show exp(n+1)           #n =1.0 t          exp(n
@show exp(p+1)           #p =  0            exp(p
@show exp(r+1)           #r =  1            exp(r
@show exp(s)+1           #s =  1 t          exp(s
@show exp(u+1)           #u =  1            exp(u
@show exp(v+1)           #v =  1 + 1 t      exp(v
@show exp(w+1)           #w =  1            exp(w
@show exp(x+1)           #x =  1            exp(x
@show exp(y+1)           #y =  1 + 1 t      exp(y 
@show exp(z+1)           #z =  1 + 1 t      exp(z
println("----------------Tutorial 5: high order Usage +1 --------------")
@show exp(nh+1)        #nh= 1.0 t           exp(nh
@show exp(ph+1)        #ph= 0               exp(ph
@show exp(rh+1)        #rh= 1               exp(rh
@show exp(sh+1)        #sh= 1 t             exp(sh
@show exp(uh+1)        #uh= 1               exp(uh
@show exp(vh+1)        #vh= 1 + 1 t         exp(vh
@show exp(wh+1)        #wh= 1               exp(wh
@show exp(xh+1)        #xh= 1 + 1 t²        exp(xh
@show exp(yh+1)        #yh= 1 + 1 t         exp(yh
@show exp(zh+1)        #zh= 1 + 1 t + 1 t²  exp(zh
=#
#= println("----------------Tutorial 6:low order sqrt--------------")
#@show sqrt(n) or sqrt(s) #illegal!!!!First non-vanishing Taylor1 coefficient must correspond to an **even power** in order to expand `sqrt` around 0
# 0.0
@show sqrt(p)           #p =  0           
 # 1.0
@show sqrt(r)           #r =  1            
@show sqrt(u)           #u =  1            
@show sqrt(w)           #w =  1           
@show sqrt(x)           #x =  1            
# 1.0 + 0.5 t
@show sqrt(v)           #v =  1 + 1 t     
@show sqrt(y)           #y =  1 + 1 t      
@show sqrt(z)  =#  
#= println("----------------Tutorial7:high order sqrt--------------")
#@show sqrt(n)# or sqrt(s) #illegal!!!!First non-vanishing Taylor1 coefficient must correspond to an **even power** in order to expand `sqrt` around 0
# 0.0
@show sqrt(ph)           #p =  0           
 # 1.0
@show sqrt(rh)           #r =  1            
@show sqrt(uh)           #u =  1            
@show sqrt(wh)           #w =  1           
#1.0 + 0.5 t² - 0.125 t⁴ 
@show sqrt(xh)           #x =  1            
# 1.0 + 0.5 t - 0.125 t² + 0.0625 t³ - 0.0390625 t⁴ 
@show sqrt(vh)           #v =  1 + 1 t     
@show sqrt(yh)           #y =  1 + 1 t      
#1.0 + 0.5 t + 0.375 t² - 0.1875 t³ + 0.0234375 t⁴ 
@show sqrt(zh)   =#


#= println("----------------Tutorial8:differentiate integrate--------------")
 @show differentiate(nh)
@show differentiate(ph)
@show differentiate(rh)
@show differentiate(sh)
@show differentiate(uh)
@show differentiate(vh)
@show differentiate(wh)
@show differentiate(xh)
@show differentiate(yh)
@show differentiate(zh) 

 @show differentiate(exp(nh))
@show differentiate(exp(ph))
@show differentiate(exp(rh))
@show differentiate(exp(sh))
@show differentiate(exp(uh))
@show differentiate(exp(vh))
@show differentiate(exp(wh))
@show differentiate(exp(xh))
@show differentiate(exp(yh))
@show differentiate(exp(zh)) 

#@show integrate((exp(nh)))
#@show integrate(differentiate(exp(nh)))
 @show differentiate(exp(ph))
@show differentiate(exp(rh))
@show differentiate(exp(sh))
@show differentiate(exp(uh))
@show differentiate(exp(vh))
@show differentiate(exp(wh))
@show differentiate(exp(xh))
@show differentiate(exp(yh))
@show differentiate(exp(zh)) 

#@show differentiate(exp(x))
#@show integrate(exp(x))=#
x = Taylor1([1,1.5,0.8],2)
y=Taylor1([0,0.0,0.0],2)
@show x
@show x+0
#= @show differentiate(x)
@show integrate(x,1)
integrate!(y,x,1)
@show y =#
#= 
 println("-----------------Tutorial 9: Type----------")
n1=Taylor1(Float64,l)
println(typeof(n1))
n2=1 +xh
println(typeof(n2)) =#
#= println("-----------------Tutorial 10: f(n1) vs f(n+1)-----------")
@show exp(zh+1)
zh1=1 +zh
@show exp(zh1) =#



#= println("-----------------Tutorial 11: evaluate-----------")

@show x(0.75)
@show y(0.75)
@show evaluate(x, 0.75)
@show evaluate(y, 0.75) =#

#= poly=Taylor1([1, 2, 3, 0.5],3)
@show poly
@show differentiate(poly)
@show integrate(poly)
@show evaluate(poly, 2.1)
@show poly(2.1) =#

#= println("-----------------Tutorial 12: any function-----------")

u33= Taylor1([1,0,1],6)
@show u33
f(t)=t*t
@show f(u33) =#

#= println("-----------------Tutorial 12: taylor and fractions-----------")
g=a-> 1/(1-a) =#

#= @show g(n)  # 1.0 + 1.0 t 
@show g(p) # 1.0
@show g(s)  # 1.0 + 1.0 t  =#
#=error
@show g(r)
@show g(u)
@show g(v)
@show g(w)
@show g(x)
@show g(y)
@show g(z)=#

#= @show g(nh)  # 1.0 + 1.0 t + 1.0 t² + 1.0 t³ + 1.0 t⁴  
@show g(ph) # 1.0
@show g(sh)  # 1.0 + 1.0 t + 1.0 t² + 1.0 t³ + 1.0 t⁴   =#




#=error
@show g(r)
@show g(u)
@show g(v)
@show g(w)
@show g(x)
@show g(y)
@show g(z)=#

