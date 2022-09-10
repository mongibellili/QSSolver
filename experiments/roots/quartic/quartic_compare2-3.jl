@inline function quadraticrealroots(b::F, c::F) where F <: AbstractFloat
    b, c = -0.5b, c
    Δ = muladd(b, b, -c) #b * b - c
    if Δ < -4eps(F) # Complex roots
        tuple()
    elseif Δ > 4eps(F) # Real roots
        q = b > 0 ? b + sqrt(Δ) : b - sqrt(Δ)
        (q, c / q)
    else # Double real root
        tuple(b)
    end
end
@inline function cubicmaxroot(b::F, c::F, d::F) where F <: AbstractFloat
    _3 = one(F) / 3
    _9 = one(F) / 9
    SQ3 = sqrt(3one(F))
    xₙ = -b * _3
    b²_9 = b * b * _9
    yₙ = muladd(muladd(-2, b²_9, c), xₙ, d) #d + xₙ * (c - 2b²_9)
    δ² = muladd(-_3, c, b²_9) #b²_9 - c *_3
    h² = 4δ² * δ² * δ²
    Δ = muladd(yₙ, yₙ, -h²) #yₙ * yₙ - h²
    if Δ > 4eps(F) # one real root and two complex roots
        p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
        q = δ² / p # cbrt(0.5 * (-yₙ - √Δ)) : cbrt(0.5 * (-yₙ + √Δ))
        xₙ + p + q
    elseif Δ < -4eps(F) # three real roots
        θ = abs(yₙ) < eps(F) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
        δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
        z₁ = 2δ * cos(θ)
        z₂ = muladd(-0.5, z₁, xₙ) #xₙ - 0.5z₁   # x₂ = xₙ + 2δ * cos(2π / 3 - θ)
        z₃ = SQ3 * δ * sin(θ) # x₃ = xₙ + 2δ * cos(4π / 3 - θ)
        x₁ = xₙ + z₁
        x₂ = z₂ + z₃
        x₃ = z₂ - z₃
        x = x₁ > x₂ ? x₁ : x₂
        x₃ > x ? x₃ : x
    else # double or triple real roots
        δ = cbrt(0.5yₙ)
        x₁ = xₙ + δ
        x₂ = xₙ - 2δ
        x₁ > x₂ ? x₁ : x₂
    end
end
function quarticsmallestpositiveroot(a::F, b::F, c::F, d::F, e::F) where F <: AbstractFloat
    _a = inv(a)
    b, c, d, e = b * _a, c * _a, d * _a, e * _a
    m = b < c ? b : c
    m = d < m ? d : m
    m = e < m ? e : m
    m > eps(F) && return typemax(F) # Cauchy bound
    xₙ = -0.25b
    b² = b * b
    p = muladd(-0.375, b², c)
    q = muladd(muladd(-0.5, b², 2c), xₙ, d)
    r = muladd(muladd(muladd(-0.1875, b², c), xₙ, d), xₙ, e)
    if abs(r) < 4eps(F)
        xₙ
    elseif abs(q) < 4eps(F)
        y² = quadraticrealroots(p, r)
        l = length(y²)
        l === 0 && return typemax(F)
        y₁ = y²[1] < -eps(F) ? typemax(F) : √y²[1]
        z₁ = xₙ + y₁
        z₂ = xₙ - y₁
        x = z₁ > -eps(F) ? z₁ : typemax(F)
        x = z₂ > -eps(F) && z₂ < x ? z₂ : x
        l === 1 && return x
        y₂ = y²[2] < -eps(F) ? typemax(F) : √y²[2]
        z₃ = xₙ + y₂
        z₄ = xₙ - y₂
        x = z₃ > -eps(F) && z₃ < x ? z₃ : x
        z₄ > -eps(F) && z₄ < x ? z₄ : x
    else
        h² = cubicmaxroot(2p, muladd(p, p, -4r), -q * q)
        h = √h²
        j = 0.5(p + h² - q / h)
        y = quadraticrealroots(h, j)
        l = length(y)
        x = if l === 0
            typemax(F)
        elseif l === 1
            z₁ = xₙ + y[1]
            z₁ > -eps(F) ? z₁ : typemax(F)
        else
            z₁ = xₙ + y[1]
            z₁ = z₁ > -eps(F) ? z₁ : typemax(F)
            z₂ = xₙ + y[2]
            z₂ > -eps(F) && z₂ < z₁ ? z₂ : z₁
        end
        y = quadraticrealroots(-h, r / j)
        l = length(y)
        if l === 0
            x
        elseif l === 1
            z₁ = xₙ + y[1]
            #z₁ > -eps(F) ? z₁ : x
            z₁ > -eps(F) && z₁ < x ? z₁ : x
        else
            z₁ = xₙ + y[1]
            #x = z₁ > -eps(F) ? z₁ : x
            x = z₁ > -eps(F) && z₁ < x ? z₁ : x
            z₂ = xₙ + y[2]
            z₂ > -eps(F) && z₂ < x ? z₂ : x
        end
    end
end
function quartic2(a::F, b::F, c::F, d::F, e::F) where F <: AbstractFloat
    _a = inv(a)
    b, c, d, e = b * _a, c * _a, d * _a, e * _a
    m = b < c ? b : c
    m = d < m ? d : m
    m = e < m ? e : m
    m > eps(F) && return typemax(F) # Cauchy bound
    xₙ = -0.25b
    b² = b * b
    p = muladd(-0.375, b², c)
    q = muladd(muladd(-0.5, b², 2c), xₙ, d)
    r = muladd(muladd(muladd(-0.1875, b², c), xₙ, d), xₙ, e)
    #x^4+px^2+qx+r=0
    if abs(r) < 4eps(F)  # r==0
        z₁=xₙ
        z₂,z₃,z₄=cubic(1.0,0.0,p,q)
    elseif abs(q) < 4eps(F)   #q==0 biquadratic 
        y² = quadraticrealroots(p, r)
        l = length(y²)
        l === 0 && return typemax(F)
        y₁ = y²[1] < -eps(F) ? typemax(F) : √y²[1]
        z₁ = xₙ + y₁
        z₂ = xₙ - y₁
       #=  x = z₁ > -eps(F) ? z₁ : typemax(F)
        x = z₂ > -eps(F) && z₂ < x ? z₂ : x =#
        if l === 1 
        z₃ = z₁
        z₄ = z₂
        else
        y₂ = y²[2] < -eps(F) ? typemax(F) : √y²[2]
        z₃ = xₙ + y₂
        z₄ = xₙ - y₂
        end
       #=  x = z₃ > -eps(F) && z₃ < x ? z₃ : x
        z₄ > -eps(F) && z₄ < x ? z₄ : x =#
    else
        h² = cubicmaxroot(2p, muladd(p, p, -4r), -q * q)
        h = √h²
        j = 0.5(p + h² - q / h)
        y = quadraticrealroots(h, j)
        l = length(y)
        #x = 
        if l === 0
            z₁,z₂=Inf,Inf#typemax(F)
        elseif l === 1
            z₁ = xₙ + y[1]
            #z₁ > -eps(F) ? z₁ : typemax(F)
            z₂=Inf
        else
            z₁ = xₙ + y[1]
           # z₁ = z₁ > -eps(F) ? z₁ : typemax(F)
            z₂ = xₙ + y[2]
           # z₂ > -eps(F) && z₂ < z₁ ? z₂ : z₁
        end
        y = quadraticrealroots(-h, r / j)
        l = length(y)
        if l === 0
            z₃,z₄=Inf,Inf#x
        elseif l === 1
            z₃ = xₙ + y[1]
            z₄=Inf
            #z₁ > -eps(F) ? z₁ : x
            #z₁ > -eps(F) && z₁ < x ? z₁ : x
        else
            z₃ = xₙ + y[1]
            #x = z₁ > -eps(F) ? z₁ : x
            #x = z₁ > -eps(F) && z₁ < x ? z₁ : x
            z₄ = xₙ + y[2]
            #z₂ > -eps(F) && z₂ < x ? z₂ : x
        end
    end
    return z₁,z₂,z₃,z₄
end
function quartic3(a::F, b::F, c::F, d::F, e::F) where F <: AbstractFloat
    _a = inv(a)
    b, c, d, e = b * _a, c * _a, d * _a, e * _a
    m = b < c ? b : c
    m = d < m ? d : m
    m = e < m ? e : m
    m > eps(F) && return typemax(F) # Cauchy bound
    xₙ = -0.25b
    b² = b * b
    p = muladd(-0.375, b², c)
    q = muladd(muladd(-0.5, b², 2c), xₙ, d)
    r = muladd(muladd(muladd(-0.1875, b², c), xₙ, d), xₙ, e)
    #x^4+px^2+qx+r=0
    if abs(r) < 4eps(F)  # r==0
        z₁=xₙ
        z₂,z₃,z₄=cubic(1.0,0.0,p,q)
    elseif abs(q) < 4eps(F)   #q==0 biquadratic 
        y² = quadraticrealroots(p, r)
        l = length(y²)
        l === 0 && return typemax(F)
        y₁ = y²[1] < -eps(F) ? typemax(F) : √y²[1]
        z₁ = xₙ + y₁
        z₂ = xₙ - y₁
        if l === 1 
        z₃ = z₁
        z₄ = z₂
        else
        y₂ = y²[2] < -eps(F) ? typemax(F) : √y²[2]
        z₃ = xₙ + y₂
        z₄ = xₙ - y₂
        end
    else
       # h² = cubicmaxroot(2p, muladd(p, p, -4r), -q * q)
       h² = cubicmaxroot(-p/2, -r,p*r/2 -q * q/8)
       # h = √h²
       h=sqrt(2h²-p)
       # j = 0.5(p + h² - q / h)
       j=-q/(2*h)+h²
        y = quadraticrealroots(h, j)
        l = length(y)
        #x = 
        if l === 0
            z₁,z₂=Inf,Inf#typemax(F)
        elseif l === 1
            z₁ = xₙ + y[1]
            z₂=Inf
        else
            z₁ = xₙ + y[1]
            z₂ = xₙ + y[2]
        end
        j=q/(2*h)+h²
        y = quadraticrealroots(-h,  j)
        l = length(y)
        if l === 0
            z₃,z₄=Inf,Inf#x
        elseif l === 1
            z₃ = xₙ + y[1]
            z₄=Inf
        else
            z₃ = xₙ + y[1]
            z₄ = xₙ + y[2]
        end
    end
    return z₁,z₂,z₃,z₄
end
function cubic(a::F, b::F, c::F, d::F) where F <: AbstractFloat
    _a = one(F) / a
    b, c, d = b * _a, c * _a, d * _a
   #=  m = b < c ? b : c
    m = d < m ? d : m
    m > eps(F) && return typemax(F) # Cauchy bound =#
    _3 = one(F) / 3
    _9 = one(F) / 9
    SQ3 = sqrt(3one(F))
    xₙ = -b * _3
    b²_9 = b * b * _9
    yₙ = muladd(muladd(-2, b²_9, c), xₙ, d)   #eq to 2R
    δ² = muladd(-_3, c, b²_9)                  #eq to Q
    h² = 4δ² * δ² * δ²
    Δ = muladd(yₙ, yₙ, -h²)
    if Δ > 4eps(F) # one real root and two complex roots
      p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
      q = δ² / p
      res1 = xₙ + p + q
      res2=Inf
      res3=Inf
      #z > -eps(F) ? z : typemax(F)
    elseif Δ < -4eps(F) # three real roots
      θ = abs(yₙ) < eps(F) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
      δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
      z₁ = 2δ * cos(θ)
      z₂ = muladd(-0.5, z₁, xₙ)
      z₃ = SQ3 * δ * sin(θ)
      res1 = xₙ + z₁
      res2 = z₂ + z₃
      res3 = z₂ - z₃
      #= x = x₁ > -eps(F) ? x₁ : typemax(F)
      x = x₂ > -eps(F) && x₂ < x ? x₂ : x
      x₃ > -eps(F) && x₃ < x ? x₃ : x =#
    else # double repeated real roots
      δ = cbrt(0.5yₙ)
      res1 = xₙ + δ
      res2 = xₙ - 2δ
      res3 =Inf
    #=  x = x₁ > -eps(F) ? x₁ : typemax(F)
      x₂ > -eps(F) && x₂ < x ? x₂ : x =#
    end
    return res1,res2,res3
end 
function quartic1(a::F, b::F, c::F, d::F, e::F) where F <: AbstractFloat
    if a == 0.0
        res1=Inf
        res2,res3,res4 = cubic(b,c,d,e)
    else
        #res = Array{Complex{Float64}}(4)
        if e == 0.0
        res1 = 0.0
        res2,res3,res4 = cubic(a,b,c,d)
        else
        a₀ = e/a
        a₁ = d/a
        a₂ = c/a
        a₃ = b/a
        y₁ = cubic(1.0,-a₂, a₁*a₃-4a₀,  4a₂*a₀-a₁^2-a₃^2*a₀)[1]
        R = sqrt(0.25a₃^2-a₂+y₁)
        if R == 0.0
            A = 0.75a₃^2-2a₂
            B = 2*sqrt(y₁^2-4a₀)
        else
            A = 0.75a₃^2-R^2-2a₂
            B = (a₃*a₂-2a₁-0.25a₃^3)/R
        end
        D = sqrt(A+B)
        E = sqrt(A-B)
        res1 = -0.25a₃+0.5R+0.5D
        res2 = res[1]-D
        res3 = -0.25a₃-0.5R+0.5E
        res4 = res[3]-E
        end
    end
    return res1,res2,res3,res4
end  # error of sqrt in this function
  
using BenchmarkTools
#pol(x)=-1.33x^4+2.0x^3+1.0x^2+6.0x+2.2
#@show pol(2.573862048253818)
#@show quartic1(-1.33,2.0,1.0,6.0,2.2)#error sqrt (-)
#= @show quartic2(-1.33,2.0,1.0,6.0,2.2) #2.573862048253818
@show quartic3(-1.33,2.0,1.0,6.0,2.2) #2.573862048253818 =#
#@show quarticsmallestpositiveroot(-1.33,2.0,1.0,6.0,2.2) #2.573862048253818

#= @btime quartic2(-1.33,2.0,1.0,6.0,2.2) #2.573862048253818
@btime quartic3(-1.33,2.0,1.0,6.0,2.2) #2.573862048253818 =#
#@btime quarticsmallestpositiveroot(-1.33,2.0,1.0,6.0,2.2) #2.573862048253818

#= c = vec(3(rand(1,3) .- 0.5))
@show c =#
#= using Random
#= r = rand(1:15)
@show r
Random.seed!(r) =#
@show rand(1,4) =#





















#testing equality of cubic equations

#= p=0.5
q=0.33
r=-3.7
pol1(x)=x^3+2p*x^2+(p*p-4r)x-q*q  #x=2y-p
pol11(x)=(2x-p)^3+2p*(2x-p)^2+(p*p-4r)*(2x-p)-q*q  #x=2y-p
pol2(x)=8*(x^3-0.5p*x^2-r*x+p*r/2-q*q/8) #x=(y+p)/2
x=2*2.2-p
@show x
@show pol1(x)
@show pol11(2.2)
@show pol2(2.2) =#