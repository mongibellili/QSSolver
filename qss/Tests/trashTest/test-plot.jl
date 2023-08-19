using Plots
using StaticArrays
#using TaylorSeries
#= function test1()
a=20;x=1;u=-20.005;Δ=x*1e-2
@show a*x+u-a*Δ
f(h)=(a*x+u)*h/(1-h*a)
h1=Δ/(a*(x+Δ)+u)
h2=-Δ/(a*(x-Δ)+u)
@show h1,h2
@show f(h1),f(h2)#= ,f(1),f(5),f(10),f(20) =#
temp=1/a+(x+u/a)/(a*Δ-a*x-u)
@show temp
@show f(temp)
end =#


#= function test2()
    a=-20;x=1;u=22.005;Δ=x*1e-2
    @show abs(a*x+u)-abs(a*Δ)
    f(h)=(a*x+u)*h/(1-h*a)
   
    h1=Δ/(a*(x+Δ)+u)
    h2=-Δ/(a*(x-Δ)+u)
    @show h1,h2
    @show f(h1),f(h2)#= ,f(1),f(5),f(10),f(20) =#
    temp=1/a+(x+u/a)/(a*Δ-a*x-u)
    @show temp,f(temp)
    @show f(1e-6+1/a)
    display(plot!(f,xlims=(-2,2),ylims=(-10.1,10.1)))
readline()
println("press keyboard")
end =#
#= function fun(h)
  aii=-20;aij=-80.0;uij=1600.0;ajj=1.24;aji=-0.01;uji=0.2;xi=1;xj=-1.0;u=22.005;Δi=xi*1e-2
  Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if abs(Δ)<1e-1
          Δ=1
        end
  return -(h*(aii*xi+aij*xj+uij)+h*h*(aij*(xi*aji+uji)-ajj*(aii*xi+uij)))/Δ

end
function test3()
    

   
  #=   h1=Δ/(a*(x+Δ)+u)
    h2=-Δ/(a*(x-Δ)+u)
    @show h1,h2
    @show f(h1),f(h2)#= ,f(1),f(5),f(10),f(20) =#
    temp=1/a+(x+u/a)/(a*Δ-a*x-u)
    @show temp,f(temp)
    @show f(1e-6+1/a) =#
    #= t=Taylor1(Float64,9)
    g=derivative(f(t))
    @show g =#
   # gg(t)=1660.0 - 66232.0*t + 1.99125648e6*t*t - 5.319904048640001e7*t^3 + 1.3324777758860803e9*t^4 - 3.20395828454252e10*t^5 + 7.4899567336463e11*t^6 - 1.7152084998061518e13*t^7 + 3.866474114759756e14*t^8
   # display(plot!(f,xlims=(-2,10),ylims=(-100.1,100.1)))
    display(plot!(fun,xlims=(-2,10),ylims=(-100.1,100.1)))
readline()
println("press keyboard")
end =#
#test3()
function PosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
  mpr=() #coef1=c, coef2=b, coef3=a
  a=coeff[3];b=coeff[2];c=coeff[1]
  if a == 0 #|| (10000 * abs(coeff[3])) < abs(coeff[2])# coef3 is the coef of t^2
      if b != 0
        if -c / b>0
         mpr = (-c / b,)
        end
      end
  else 
     #double disc;
     Δ = 1.0 - 4c*a / (b*b)
        if Δ >0
          q = -0.5*(1.0+sign(b)*sqrt(Δ))*b
          r1 = q / a
         
          r2=c / q
         
       
        @show r1,r2
        if ((r1 > 0) && (r2 >0)) 
          mpr = (r1,r2)
        elseif ((r1 > 0) && (r2 <=0)) 
          mpr = (r1,)
        elseif ((r1 < 0) && (r2 >0)) 
          mpr = (r2,)
        end
      end
  end
  return mpr
end
#= function PosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
  mpr=() #coef1=c, coef2=b, coef3=a
  a=coeff[3];b=coeff[2];c=coeff[1]
  @show a,b,c
  if a == 0 #|| (10000 * abs(coeff[3])) < abs(coeff[2])# coef3 is the coef of t^2
      if b != 0
        if -c / b>0
         mpr = (-c / b,)
        end
      end
  else 
     #double disc;
     Δ = b*b-4*a*c
        if Δ >0
          @show sqrt(Δ),b
          r1 = (-b+sqrt(Δ))/(2*a)
         
          r2=(-b-sqrt(Δ))/(2*a)
         
       
        @show r1,r2
        if ((r1 > 0) && (r2 >0)) 
          mpr = (r1,r2)
        elseif ((r1 > 0) && (r2 <=0)) 
          mpr = (r1,)
        elseif ((r1 < 0) && (r2 >0)) 
          mpr = (r2,)
        end
      end
  end
  return mpr
end =#
function plotfj(a,b,c)
  fh(x)=c +b*x+ a*x*x
  #= @show fh(1.6*1e-12)
  @show fh(1.3378773055158425e-12)
  @show fh(50876.52454563779)
 # fh2(x)=-1.0e-6-747454.6988035275x+ 14.651504647354596*x*x
  coefj=@SVector [1.0e-6, -747452.6968034876, 14.691504647354595]
  posSolPlusj= PosRoot(coefj, Val(2))
  @show posSolPlusj =#
  #display(plot!(fh,xlims=(-1e-11,1e-11),ylims=(-0.000001,0.000001)))
  #display(plot!(fh))
  display(plot!(fh,xlims=(0,0.01),ylims=(-1.0,1.0)))
  readline()
  println("press keyboard")
end
#= a,b,c=-185269.7656774946, -102.3518745108343, -0.014735813027671029
plotfj(a,b,c) =#



#= resppi = pointer(Vector{Float64}(undef, 2))
@show ccc =#