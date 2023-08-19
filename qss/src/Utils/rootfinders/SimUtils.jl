using StaticArrays

function minPosRoot(coeff::SVector{2,Float64}, ::Val{1}) # coming from val(1) means coef has x and derx only...size 2
    mpr=-1
        if coeff[2] == 0 
            mpr = Inf
        else 
            mpr = -coeff[1] / coeff[2];
        end
        if mpr < 0
            mpr = Inf
        end
       # println("mpr inside minPosRoot in utils= ",mpr)
    return mpr
end
function minPosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
    mpr=-1 #coef1=c, coef2=b, coef3=a
    if coeff[3] == 0 #|| (10000 * abs(coeff[3])) < abs(coeff[2])# coef3 is the coef of t^2
        if coeff[2] == 0
          mpr = Inf
        else 
          mpr = -coeff[1] / coeff[2]
        end
        if mpr < 0
          mpr = Inf
        end
    else 
       #double disc;
        disc = coeff[2] * coeff[2] - 4 * coeff[3] * coeff[1]#b^2-4ac
        if disc < 0 # no real roots
          mpr = Inf
        else 
          #double sd, r1;
          sd = sqrt(disc);
          r1 = (-coeff[2] + sd) / (2 * coeff[3]);
          if r1 > 0 
            mpr = r1;
          else 
            mpr = Inf;
          end
          r1 = (-coeff[2] - sd) / (2 * coeff[3]);
          if ((r1 > 0) && (r1 < mpr)) 
            mpr = r1;
          end
        end
        
    end
    return mpr
end

#= function minPosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
  mpr=-1 #coef1=c, coef2=b, coef3=a
  a=coeff[3];b=coeff[2];c=coeff[1]
  if a == 0 #|| (10000 * abs(coeff[3])) < abs(coeff[2])# coef3 is the coef of t^2
      if b == 0
        mpr = Inf
      else 
        mpr = -c / b
      end
      if mpr < 0
        mpr = Inf
      end
  else 
     #double disc;
     Δ = 1.0 - 4c*a / (b*b)
     if Δ < 0.0
      mpr=Inf
     else
       q = -0.5*(1.0+sign(b)*sqrt(Δ))*b
       mpr = q / a
       if mpr<=0
         mpr=Inf
       end
       mpr2=c / q
       if  mpr2>0 && mpr2<mpr
        mpr=mpr2
       end
     end
      
  end
  return mpr
end
 =#

#= @inline function minPosRoot(coeffs::SVector{4,Float64}, ::Val{3})#where F <: AbstractFloat
  if coeffs[4] == 0.0
    coeffs2=@SVector[coeffs[1],coeffs[2],coeffs[3]]
    return minPosRoot(coeffs2, Val(2))
  end
  _a = 1.0 / coeffs[4]
  b, c, d = coeffs[3] * _a, coeffs[2] * _a, coeffs[1] * _a
  m = b < c ? b : c
  m = d < m ? d : m
  m > eps(Float64) && return Inf#typemax(Float64) # Cauchy bound
  _3 = 1.0 / 3
  _9 = 1.0 / 9
  SQ3 = sqrt(3.0)
  xₙ = -b * _3
  b²_9 = b * b * _9
  yₙ = muladd(muladd(-2, b²_9, c), xₙ, d)   #eq to 2R
  δ² = muladd(-_3, c, b²_9)                  #eq to Q
  h² = 4δ² * δ² * δ²
  Δ = muladd(yₙ, yₙ, -h²)
  if Δ > 4eps(Float64) # one real root and two complex roots
  p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
  q = δ² / p
  z = xₙ + p + q
  z > -eps(Float64) ? z : Inf#typemax(Float64)
  elseif Δ < -4eps(Float64) # three real roots
  θ = abs(yₙ) < eps(Float64) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
  δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
  z₁ = 2δ * cos(θ)
  z₂ = muladd(-0.5, z₁, xₙ)
  z₃ = SQ3 * δ * sin(θ)
  x₁ = xₙ + z₁
  x₂ = z₂ + z₃
  x₃ = z₂ - z₃
  x = x₁ > -eps(Float64) ? x₁ : Inf# typemax(F)
  x = x₂ > -eps(Float64) && x₂ < x ? x₂ : x
  x₃ > -eps(Float64) && x₃ < x ? x₃ : x
  else # double or triple real roots
  δ = cbrt(0.5yₙ)
  x₁ = xₙ + δ
  x₂ = xₙ - 2δ
  x = x₁ > -eps(Float64) ? x₁ : Inf#typemax(F)
  x₂ > -eps(Float64) && x₂ < x ? x₂ : x
  end
end =#

@inline function minPosRoot(coeffs::SVector{4,Float64}, ::Val{3})#where F <: AbstractFloat
  if coeffs[4] == 0.0
    coeffs2=@SVector[coeffs[1],coeffs[2],coeffs[3]]
    return minPosRoot(coeffs2, Val(2))
  end
  _a = 1.0 / coeffs[4]
  b, c, d = coeffs[3] * _a, coeffs[2] * _a, coeffs[1] * _a
  m = b < c ? b : c
  m = d < m ? d : m
  m > 0.0 && return Inf#typemax(Float64) # Cauchy bound
  _3 = 1.0 / 3
  _9 = 1.0 / 9
  SQ3 = sqrt(3.0)
  xₙ = -b * _3
  b²_9 = b * b * _9
  yₙ = muladd(muladd(-2, b²_9, c), xₙ, d)   #eq to 2R
  δ² = muladd(-_3, c, b²_9)                  #eq to Q
  h² = 4δ² * δ² * δ²
  Δ = muladd(yₙ, yₙ, -h²)
  if Δ > 0.0 # one real root and two complex roots
  p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
  q = δ² / p
  z = xₙ + p + q
  z > 0.0 ? z : Inf#typemax(Float64)
  elseif Δ < 0.0 # three real roots
  θ = abs(yₙ) < 0.0 ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
  δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
  z₁ = 2δ * cos(θ)
  z₂ = muladd(-0.5, z₁, xₙ)
  z₃ = SQ3 * δ * sin(θ)
  x₁ = xₙ + z₁
  x₂ = z₂ + z₃
  x₃ = z₂ - z₃
  x = x₁ > 0.0 ? x₁ : Inf# typemax(F)
  x = x₂ > 0.0 && x₂ < x ? x₂ : x
  x₃ > 0.0 && x₃ < x ? x₃ : x
  else # double or triple real roots
  δ = cbrt(0.5yₙ)
  x₁ = xₙ + δ
  x₂ = xₙ - 2δ
  x = x₁ > 0.0 ? x₁ : Inf#typemax(F)
  x₂ > 0.0 && x₂ < x ? x₂ : x
  end
end
###########later optimize

#= @inline function minPosRoot(coeffs::NTuple{3,Float64}, ::Val{2}) 


	_a = 1.0 / coeffs[1]
	b, c = -0.5coeffs[2] * _a, coeffs[3] * _a
	Δ = muladd(b, b, -c) # b * b - c
	if Δ < -4eps() # Complex roots
		Inf
	elseif Δ > 4eps() # Real roots
		if b > eps()
			c > eps() ? c / (b + sqrt(Δ)) : b + sqrt(Δ)
		elseif b < -eps()
			c < eps() ? c / (b - sqrt(Δ)) : Inf
		else
			sqrt(-c)
		end
	else # Double real root
		b > -eps() ? b : Inf
	end
end =#

#coef=@SVector [-5.144241008311048e7, -8938.3815305787700, -0.2906652780025, 1e-6]
#= coef=@SVector [1e-6, -0.2906652780025, -8938.3815305787700, -5.144241008311048e7]
x=minPosRoot(coef, Val(3))
@show x =#
#coef=@SVector [-5.144241008311048e7, -8938.3815305787700, -0.2906652780025, 1e-6]
#= coef=@SVector [1e-6, -0.040323840000000166, -3914.116824448214, 1.021243315770821e8]
x=minPosRoot(coef, Val(3))
@show x =#


#= function PosRoot(coeff::SVector{3,Float64}, ::Val{2}) # credit goes to github.com/CIFASIS/qss-solver
  mpr=() #coef1=c, coef2=b, coef3=a
  if coeff[3] == 0 #|| (10000 * abs(coeff[3])) < abs(coeff[2])# coef3 is the coef of t^2
      if coeff[2] != 0
        if -coeff[1] / coeff[2]>0
         mpr = (-coeff[1] / coeff[2],)
        end
      end
  else 
     #double disc;
      disc = coeff[2] * coeff[2] - 4 * coeff[3] * coeff[1]#b^2-4ac
      if disc >0
        #double sd, r1;
        sd = sqrt(disc);
        r1 = (-coeff[2] + sd) / (2 * coeff[3]);
        r2 = (-coeff[2] - sd) / (2 * coeff[3]);
        if ((r1 > 0) && (r2 >0)) 
          mpr = (r1,r2)
        elseif ((r1 > 0) && (r2 <0)) 
          mpr = (r1,)
        elseif ((r1 < 0) && (r2 >0)) 
          mpr = (r2,)
        end
      end
  end
  return mpr
end =#

function PosRoot(coeff::SVector{3,Float64}, ::Val{2}) # 
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
         
       
        #@show r1,r2
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

#= function sort3(tup)
  h1=min(tup[1],tup[2],tup[3])
  h3=max(tup[1],tup[2],tup[3])
  if tup[1]!=h1 && tup[1]!=h3
    h2=tup[1]
  elseif tup[2]!=h1 && tup[2]!=h3
    h2=tup[2]
  else
    h2=tup[3]
  end
  res=(h1,h2,h3)
end

constructIntrval(tup::Tuple{})=((0,Inf),)
constructIntrval(tup::Tuple{Float64})=((0,tup[1]),)
function constructIntrval(tup::Tuple{Float64,Float64})
    h1=min(tup[1],tup[2])
    h4=max(tup[1],tup[2])
    res=((0,h1),(h4,Inf))
end
function constructIntrval(tup::NTuple{3, Float64})
    h1,h2,h3=sort3(tup)
    res=((0,h1),(h2,h3))
end
function constructIntrval(tup::NTuple{4, Float64})
    h1=min(tup[1],tup[2],tup[3],tup[4])
    if tup[1]==h1
      h2,h3,h4=sort3((tup[2],tup[3],tup[4]))
    elseif tup[2]==h1
      h2,h3,h4=sort3((tup[1],tup[3],tup[4]))
    elseif tup[3]==h1
      h2,h3,h4=sort3((tup[2],tup[1],tup[4]))
    else
      h2,h3,h4=sort3((tup[2],tup[3],tup[1]))
    end
    #if c/β << delta
      #solve for 
    res=((0,h1),(h2,h3),(h4,Inf))
 
end =#

#= 
#merge is empty
constructIntrval(tup::Tuple{},tup2::Tuple{})=((0,Inf),)

#merge has 1 element
constructIntrval(tup::Tuple{Float64},tup2::Tuple{})=((0,tup[1]),)
constructIntrval(tup::Tuple{},tup2::Tuple{Float64})=constructIntrval(tup2,tup)

#merge has 2 elements
function constructIntrval(tup::Tuple{Float64},tup2::Tuple{Float64})
  h1=min(tup[1],tup2[1])
  h4=max(tup[1],tup2[1])
  res=((0,h1),(h4,Inf))
end
function constructIntrval(tup::Tuple{Float64,Float64},tup2::Tuple{})
    h1=min(tup[1],tup[2])
    h4=max(tup[1],tup[2])
    res=((0,h1),(h4,Inf))
end
constructIntrval(tup2::Tuple{},tup::Tuple{Float64,Float64})=constructIntrval(tup,tup2)

#merge has 3 elements
function constructIntrval(tup::NTuple{2, Float64},tup2::Tuple{Float64})
  t1=tup[1];t2=tup[2];t3=tup2[1]
  t12=min(t1,t2);t21=max(t1,t2)
  if t3<t12
   res=((0,t3),(t12,t21))
  else
    if t3<t21
      res=((0,t12),(t3,t21))
    else
      res=((0,t12),(t21,t3))
    end
  end
  res
end
constructIntrval(tup2::Tuple{Float64},tup::NTuple{2, Float64})=constructIntrval(tup,tup2)

#merge has 4 elements
function constructIntrval(tup::NTuple{2, Float64},tup2::NTuple{2, Float64})
  t1=tup[1];t2=tup[2];t3=tup2[1];t4=tup2[2]
    h1=min(t1,t2,t3,t4);h4=max(t1,t2,t3,t4)
    if t1==h1
      if t2==h4
        res=((0,t1),(min(t3,t4),max(t3,t4)),(t2,Inf))
      elseif t3==h4
        res=((0,t1),(min(t2,t4),max(t2,t4)),(t3,Inf))
      else
        res=((0,t1),(min(t2,t3),max(t2,t3)),(t4,Inf))
      end

    elseif t2==h1
      if t1==h4
        res=((0,t2),(min(t3,t4),max(t3,t4)),(t1,Inf))
      elseif t3==h4
        res=((0,t2),(min(t1,t4),max(t1,t4)),(t3,Inf))
      else
        res=((0,t2),(min(t1,t3),max(t1,t3)),(t4,Inf))
      end
     
    elseif t3==h1
      if t1==h4
        res=((0,t3),(min(t2,t4),max(t2,t4)),(t1,Inf))
      elseif t2==h4
        res=((0,t3),(min(t1,t4),max(t1,t4)),(t2,Inf))
      else
        res=((0,t3),(min(t1,t2),max(t1,t2)),(t4,Inf))
      end
      
    else
      if t1==h4
        res=((0,t4),(min(t2,t3),max(t2,t3)),(t1,Inf))
      elseif t2==h4
        res=((0,t4),(min(t1,t3),max(t1,t3)),(t2,Inf))
      else
        res=((0,t4),(min(t1,t2),max(t1,t2)),(t3,Inf))
      end
    end
   
    res
 
end =#



#= constructIntrval0(tup::Tuple{},tup2::Tuple{})=((0,Inf),)

#merge has 1 element
function constructIntrval1(tup::Tuple{Float64},tup2::NTuple{0, Float64})
    res=((0,tup[1]),)
end
function constructIntrval1(tup::NTuple{0, Float64},tup2::NTuple{1, Float64})
    #constructIntrval(tup2,tup)
    res=((0,tup2[1]),)
end

#merge has 2 elements
function constructIntrval2(tup::Tuple{Float64},tup2::Tuple{Float64})
  h1=min(tup[1],tup2[1])
  h4=max(tup[1],tup2[1])
  res=((0.0,h1),(h4,Inf))
end

function constructIntrval2(tup::NTuple{2, Float64},tup2::NTuple{0, Float64})
    h1=min(tup[1],tup[2])
    h4=max(tup[1],tup[2])
    res=((0.0,h1),(h4,Inf))
end 

function constructIntrval2(tup2::NTuple{0, Float64},tup::NTuple{2, Float64})
    h1=min(tup[1],tup[2])
    h4=max(tup[1],tup[2])
    res=((0.0,h1),(h4,Inf))
end 
#constructIntrval(tup2::Tuple{},tup::NTuple{2, Float64})=constructIntrval(tup,tup2)

#merge has 3 elements
function constructIntrval3(tup::NTuple{2, Float64},tup2::NTuple{1,Float64})
  t1=tup[1];t2=tup[2];t3=tup2[1]
  #t12=min(t1,t2);t21=max(t1,t2)
  if t1<t2
    t12=t1;t21=t2
  else
    t12=t2;t21=t1
  end

  if t3<t12
   res=((0.0,t3),(t12,t21))
  else
    if t3<t21
      res=((0.0,t12),(t3,t21))
    else
      res=((0.0,t12),(t21,t3))
    end
  end
  res
end
function constructIntrval3(tup2::Tuple{Float64},tup::NTuple{2, Float64})
   # constructIntrval(tup,tup2)
   t1=tup[1];t2=tup[2];t3=tup2[1]
  #t12=min(t1,t2);t21=max(t1,t2)
  if t1<t2
    t12=t1;t21=t2
  else
    t12=t2;t21=t1
  end

  if t3<t12
   res=((0.0,t3),(t12,t21))
  else
    if t3<t21
      res=((0.0,t12),(t3,t21))
    else
      res=((0.0,t12),(t21,t3))
    end
  end
  res
end

#merge has 4 elements
function constructIntrval4(tup::NTuple{2, Float64},tup2::NTuple{2, Float64})
  t1=tup[1];t2=tup[2];t3=tup2[1];t4=tup2[2]
    h1=min(t1,t2,t3,t4);h4=max(t1,t2,t3,t4)
    if t1==h1
      if t2==h4
        res=((0,t1),(min(t3,t4),max(t3,t4)),(t2,Inf))
      elseif t3==h4
        res=((0,t1),(min(t2,t4),max(t2,t4)),(t3,Inf))
      else
        res=((0,t1),(min(t2,t3),max(t2,t3)),(t4,Inf))
      end

    elseif t2==h1
      if t1==h4
        res=((0,t2),(min(t3,t4),max(t3,t4)),(t1,Inf))
      elseif t3==h4
        res=((0,t2),(min(t1,t4),max(t1,t4)),(t3,Inf))
      else
        res=((0,t2),(min(t1,t3),max(t1,t3)),(t4,Inf))
      end
     
    elseif t3==h1
      if t1==h4
        res=((0,t3),(min(t2,t4),max(t2,t4)),(t1,Inf))
      elseif t2==h4
        res=((0,t3),(min(t1,t4),max(t1,t4)),(t2,Inf))
      else
        res=((0,t3),(min(t1,t2),max(t1,t2)),(t4,Inf))
      end
      
    else
      if t1==h4
        res=((0,t4),(min(t2,t3),max(t2,t3)),(t1,Inf))
      elseif t2==h4
        res=((0,t4),(min(t1,t3),max(t1,t3)),(t2,Inf))
      else
        res=((0,t4),(min(t1,t2),max(t1,t2)),(t3,Inf))
      end
    end
   
    res
 
end



#= 
vec1=(23.3,)
vec2=() =#

#display(@code_warntype constructIntrval(vec1,vec2))
#= @show constructIntrval0(vec1,vec2)
@btime constructIntrval0(vec1,vec2) =#
#@show ress

#= oi=NTuple{3, Float64}((2.0,3.2,2.2))
display(typeof(oi)) =#

#println(methods(constructIntrval))

function constructIntrval(vec1::NTuple{M, Float64},vec2::NTuple{N, Float64})where{M,N}
    #= vec1=(1.2,3.3)
    vec2=() =#
    lv1=length(vec1);lv2=length(vec2)
    lv=length(vec1)+length(vec2)
    ress=nothing
    if lv==0
        ress=constructIntrval0(vec1,vec2)
    elseif lv==1
      ress=constructIntrval1(vec1,vec2)
    elseif lv==2
      #= @timeit "2vers" =#  ress=constructIntrval2(vec1,vec2)
    elseif lv==3
        ress=constructIntrval3(vec1,vec2)
    else
        ress=constructIntrval4(vec1,vec2)
    end
    ress
end =#




function constructIntrval2(cache::Vector{Vector{Float64}},t1::Float64,t2::Float64)
  h1=min(t1,t2)
  h4=max(t1,t2)
 # res=((0.0,h1),(h4,Inf))
  cache[1][1]=0.0;cache[1][2]=h1;
  cache[2][1]=h4;cache[2][2]=Inf;
  return nothing
end 

#constructIntrval(tup2::Tuple{},tup::NTuple{2, Float64})=constructIntrval(tup,tup2)

#merge has 3 elements
function constructIntrval3(cache::Vector{Vector{Float64}},t1::Float64,t2::Float64,t3::Float64)
#t1=tup[1];t2=tup[2];t3=tup2[1]
t12=min(t1,t2);t21=max(t1,t2)
#=  if t1<t2
  t12=t1;t21=t2
else
  t12=t2;t21=t1
end =#

if t3<t12
# res=((0.0,t3),(t12,t21))
 cache[1][1]=0.0;cache[1][2]=t3;
 cache[2][1]=t12;cache[2][2]=t21;
else
  if t3<t21
   # res=((0.0,t12),(t3,t21))
    cache[1][1]=0.0;cache[1][2]=t12;
    cache[2][1]=t3;cache[2][2]=t21;
  else
    #res=((0.0,t12),(t21,t3))
    cache[1][1]=0.0;cache[1][2]=t12;
    cache[2][1]=t21;cache[2][2]=t3;
  end
end
return nothing
end


#merge has 4 elements
function constructIntrval4(cache::Vector{Vector{Float64}},t1::Float64,t2::Float64,t3::Float64,t4::Float64,)
# t1=tup[1];t2=tup[2];t3=tup2[1];t4=tup2[2]
  h1=min(t1,t2,t3,t4);h4=max(t1,t2,t3,t4)
  if t1==h1
    if t2==h4
     # res=((0,t1),(min(t3,t4),max(t3,t4)),(t2,Inf))
      cache[1][1]=0.0;cache[1][2]=t1;
      cache[2][1]=min(t3,t4);cache[2][2]=max(t3,t4);
      cache[3][1]=t2;cache[3][2]=Inf;
    elseif t3==h4
      #res=((0,t1),(min(t2,t4),max(t2,t4)),(t3,Inf))
      cache[1][1]=0.0;cache[1][2]=t1;
      cache[2][1]=min(t2,t4);cache[2][2]=max(t2,t4);
      cache[3][1]=t3;cache[3][2]=Inf;
    else
     # res=((0,t1),(min(t2,t3),max(t2,t3)),(t4,Inf))
      cache[1][1]=0.0;cache[1][2]=t1;
      cache[2][1]=min(t3,t2);cache[2][2]=max(t3,t2);
      cache[3][1]=t4;cache[3][2]=Inf;
    end

  elseif t2==h1
    if t1==h4
      #res=((0,t2),(min(t3,t4),max(t3,t4)),(t1,Inf))
      cache[1][1]=0.0;cache[1][2]=t2;
      cache[2][1]=min(t3,t4);cache[2][2]=max(t3,t4);
      cache[3][1]=t1;cache[3][2]=Inf;
    elseif t3==h4
      #res=((0,t2),(min(t1,t4),max(t1,t4)),(t3,Inf))
      cache[1][1]=0.0;cache[1][2]=t2;
      cache[2][1]=min(t1,t4);cache[2][2]=max(t1,t4);
      cache[3][1]=t3;cache[3][2]=Inf;
    else
     # res=((0,t2),(min(t1,t3),max(t1,t3)),(t4,Inf))
      cache[1][1]=0.0;cache[1][2]=t2;
      cache[2][1]=min(t3,t1);cache[2][2]=max(t3,t1);
      cache[3][1]=t4;cache[3][2]=Inf;
    end
   
  elseif t3==h1
    if t1==h4
      #res=((0,t3),(min(t2,t4),max(t2,t4)),(t1,Inf))
      cache[1][1]=0.0;cache[1][2]=t3;
      cache[2][1]=min(t2,t4);cache[2][2]=max(t2,t4);
      cache[3][1]=t1;cache[3][2]=Inf;
    elseif t2==h4
      #res=((0,t3),(min(t1,t4),max(t1,t4)),(t2,Inf))
      cache[1][1]=0.0;cache[1][2]=t3;
      cache[2][1]=min(t1,t4);cache[2][2]=max(t1,t4);
      cache[3][1]=t2;cache[3][2]=Inf;
    else
      #res=((0,t3),(min(t1,t2),max(t1,t2)),(t4,Inf))
      cache[1][1]=0.0;cache[1][2]=t3;
      cache[2][1]=min(t1,t2);cache[2][2]=max(t1,t2);
      cache[3][1]=t4;cache[3][2]=Inf;
    end
    
  else
    if t1==h4
      #res=((0,t4),(min(t2,t3),max(t2,t3)),(t1,Inf))
      cache[1][1]=0.0;cache[1][2]=t4;
      cache[2][1]=min(t3,t2);cache[2][2]=max(t3,t2);
      cache[3][1]=t1;cache[3][2]=Inf;
    elseif t2==h4
     # res=((0,t4),(min(t1,t3),max(t1,t3)),(t2,Inf))
      cache[1][1]=0.0;cache[1][2]=t4;
      cache[2][1]=min(t3,t1);cache[2][2]=max(t3,t1);
      cache[3][1]=t2;cache[3][2]=Inf;
    else
      #res=((0,t4),(min(t1,t2),max(t1,t2)),(t3,Inf))
      cache[1][1]=0.0;cache[1][2]=t4;
      cache[2][1]=min(t1,t2);cache[2][2]=max(t1,t2);
      cache[3][1]=t3;cache[3][2]=Inf;
    end
  end
 
  return nothing

end



#= 
@inline function constructIntrval(cache::Vector{Vector{Float64}},res1::Float64,res2::Float64,res3::Float64,res4::Float64)
  #= vec1=(1.2,3.3)
  vec2=() =#
  posSol=filter((x) -> x >0.0 , (res1,res2,res3,res4))

  lv1=length(vec1);lv2=length(vec2)
  lv=length(vec1)+length(vec2)
 
  if lv==0
      #res=((0.0,Inf),)
      cache[1][1]=0.0;cache[1][2]=Inf;
  elseif lv==1
      if lv1==1
          #res=((0.0,vec1[1]),)
          cache[1][1]=0.0;cache[1][2]=vec1[1];
      else
         # res=((0.0,vec2[2]),)
          cache[1][1]=0.0;cache[1][2]=vec2[1];
      end
  elseif lv==2
      if lv1==1
           constructIntrval2(cache,vec1[1],vec2[1])
      elseif lv1==2
         constructIntrval2(cache,vec1[1],vec1[2])
      else
         constructIntrval2(cache,vec2[1],vec2[2])
      end
  elseif lv==3
      if lv1==1
         constructIntrval3(cache,vec1[1],vec2[1],vec2[2])
       else
           constructIntrval3(cache,vec1[1],vec1[2],vec2[1])
      
       end
  else
     constructIntrval4(cache,vec1[1],vec1[2],vec2[1],vec2[2])
  end
  return nothing
end =#


#= @inline function constructIntrval(cache::Vector{Vector{Float64}},res1::Float64,res2::Float64,res3::Float64,res4::Float64)
  #= vec1=(1.2,3.3)
  vec2=() =#
  posSol=filter((x) -> x >0.0 , (res1,res2,res3,res4))

  lv=length(posSol)
 
  if lv==0
      #res=((0.0,Inf),)
      cache[1][1]=0.0;cache[1][2]=Inf;
  elseif lv==1
          cache[1][1]=0.0;cache[1][2]=posSol[1];
      
  elseif lv==2
     
           constructIntrval2(cache,posSol[1],posSol[2])
    
  elseif lv==3
     
         constructIntrval3(cache,posSol[1],posSol[2],posSol[3])
      
  else
     constructIntrval4(cache,posSol[1],posSol[2],posSol[3],posSol[4])
  end
  return nothing
end =#


@inline function constructIntrval(cache::Vector{Vector{Float64}},res1::Float64,res2::Float64,res3::Float64,res4::Float64)
  #= vec1=(1.2,3.3)
  vec2=() =#
  if res1<0 && res2<0 && res3<0 && res4<0
    #res=((0.0,Inf),)
    cache[1][1]=0.0;cache[1][2]=Inf;
  elseif res1>0 && res2<0 && res3<0 && res4<0
        cache[1][1]=0.0;cache[1][2]=res1;
  elseif res2>0 && res1<0 && res3<0 && res4<0
        cache[1][1]=0.0;cache[1][2]=res2;
  elseif res3>0 && res2<0 && res1<0 && res4<0
        cache[1][1]=0.0;cache[1][2]=res3;
  elseif res4>0 && res2<0 && res3<0 && res1<0
        cache[1][1]=0.0;cache[1][2]=res4;

  elseif res1>0 && res2>0 && res3<0 && res4<0
      constructIntrval2(cache,res1,res2)
  elseif res1>0 && res3>0 && res2<0 && res4<0
      constructIntrval2(cache,res1,res3)
  elseif res1>0 && res4>0 && res2<0 && res4<0
      constructIntrval2(cache,res1,res4)
  elseif res2>0 && res3>0 && res1<0 && res4<0
      constructIntrval2(cache,res2,res3)
  elseif res2>0 && res4>0 && res1<0 && res3<0
      constructIntrval2(cache,res2,res4)
  elseif res3>0 && res4>0 && res1<0 && res2<0
      constructIntrval2(cache,res3,res4)
      

   elseif res1>0 && res2>0 && res3>0 && res4<0
      constructIntrval3(cache,res1,res2,res3)
     
    elseif res1>0 && res2>0 && res4>0 && res3<0
      constructIntrval3(cache,res1,res2,res4)
    elseif res1>0 && res3>0 && res4>0 && res2<0
      constructIntrval3(cache,res1,res3,res4)
    elseif res1<0 && res3>0 && res4>0 && res2>0
      constructIntrval3(cache,res2,res3,res4)    
      
  else
     constructIntrval4(cache,res1,res2,res3,res4)
  end
  return nothing
end