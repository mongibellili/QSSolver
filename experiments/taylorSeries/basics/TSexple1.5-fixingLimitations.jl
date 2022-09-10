using TaylorSeries
#use_show_default(true)

a= Taylor1(3)
b= Taylor1([1e-120,1,0],3)
c= Taylor1([0,0,1],6)
#= @show a
@show sqrt(a+1e-120)
@show b
@show sqrt(b)
@show c
@show sqrt(c) =#
@show sqrt(a+1)



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
