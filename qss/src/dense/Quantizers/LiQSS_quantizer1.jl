function prepareAii(i::Int,a)
    return a[i][i]
end
@inline function integrateOldx(::Val{1}, x::Taylor0,olddx::MVector{O,Float64},elapsed::Float64) where{O}
    olddx[1]=x.coeffs[2]
end
  
@inline function integrateOldx(::Val{2}, x::Taylor0,olddx::MVector{2,Float64},elapsed::Float64) 
olddx[1]=x.coeffs[2]+elapsed*x.coeffs[3]*2
end
@inline function integrateOldx(::Val{3}, x::Taylor0,olddx::MVector{3,Float64},elapsed::Float64) 
olddx[1]=x.coeffs[2]+elapsed*x.coeffs[3]*2+elapsed*elapsed*x.coeffs[4]*3
end
function prepareAii(i::Int,qv::Vector{Taylor0}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64},simt::Float64)
    cacheA[1]=0.0;exactA(qv,d,cacheA, i, i, simt);a=cacheA[1]
    return a
end
function updateQ(::Val{1},opt::Options{1,DU}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64},a::Float64, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})where{DU}
   
    #cacheA[1]=0.0;exactA(qv,d,cacheA, i, i, simt);a=cacheA[1]
    q=qv[i][0];x=xv[i][0];x1=xv[i][1];
    qaux[i][1]=q
    u=x1-a*q
    dxaux[i][1]=x1
    h=0.0
    Δ=quantum[i]

    u=x1-a*q
    multiplier=opt.multiplier
    qplus=x+sign(x1)*Δ
    dxplus=a*(qplus)+u
 
    if x1*dxplus>0.0
        q=qplus
    else
        if a!=0.0
            q=-u/a
           
        else
            q=x
        end
    end
    if (q-x)>1*Δ
        q=x+1*Δ
   
    elseif (q-x)<-1*Δ
        q=x-1*Δ
  
    end
    if x1!=0.0
        if q!=x
            nextStateTime[i]=simt+abs((q-x)/x1)
        else
            nextStateTime[i]=Inf
        end
    else
        nextStateTime[i]=Inf
    end
  
  
 
    qv[i][0]=q
    return nothing
end

function updateQ(::Val{1},opt::Options{2,DU}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64},a::Float64, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})where{DU}
     
    q = qv[i][0]
    x = xv[i][0]
    x1 = xv[i][1]
    qaux[i][1] = q
    u = x1 - a * q
    dxaux[i][1] = x1
    dx=x1
    h=0.0
    Δ=quantum[i]
    multiplier=opt.multiplier
   if a !=0.0
       if dx==0.0
           dx=u+(q)*a
           if dx==0.0
               dx=1e-26
           end
       end
   #for order1 finding h is easy but for higher orders iterations are cheaper than finding exact h using a quadratic,cubic...
   #exacte for order1: h=-2Î”/(u+xa-2aÎ”) or h=2Î”/(u+xa+2aÎ”)
       h = ft-simt
       q = (x + h * u) /(1 - h * a)
      if (abs(q - x) > multiplier* Δ) 
         h = (abs( Δ / dx));
         q= (x + h * u) /(1 - h * a)
       end
       while (abs(q - x) > multiplier* Δ) 
         #h = h * 0.98*(quantum[i] / abs(q - x));
         h = h *sqrt(Δ / abs(q - x))
         q= (x + h * u) /(1 - h * a)
       end 
      
   else
       dx=u
       if dx>0.0
           q=x+Δ# 
       else
           q=x-Δ
       end
       if dx!=0
       h=(abs(Δ/dx))
       else
           h=Inf
       end
   end
   qv[i][0]=q
  # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
  nextStateTime[i]=simt+h
  if (q-x)>multiplier*Δ
    q=x+multiplier*Δ
  elseif (q-x)<-multiplier*Δ
    q=x-multiplier*Δ
    
end
   return nothing
end




function updateQ(::Val{1},opt::Options{3,DU}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64},a::Float64, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})where{DU}
    
    q = qv[i][0]
    x = xv[i][0]
    x1 = xv[i][1]
    qaux[i][1] = q
    u = x1 - a * q
    dxaux[i][1] = x1
    h = 0.0
    Δ = quantum[i]
    multiplier=opt.multiplier
    if a != 0.0
        α = -(a * x + u) / a
        h1denom = a * (x + Δ) + u
        h2denom = a * (x - Δ) + u
        h1 = Δ / h1denom
        h2 = -Δ / h2denom
        if a < 0
            if α > Δ
                h = h1
                q = x + Δ
            elseif α < -Δ
                h = h2
                q = x - Δ
            else
                h = Inf
                q = -u / a
            end
        else
            if α > Δ
                h = h2
                q = x - Δ
            elseif α < -Δ
                h = h1
                q = x + Δ
            else
                if a * x + u > 0
                    q = x - Δ
                    h = h2
                else
                    q = x + Δ
                    h = h1
                end
            end
        end
    else
        if x1 > 0.0
            q = x + Δ
        else
            q = x - Δ
        end
        if x1 != 0
            h = (abs(Δ / x1))
        else
            h = Inf
        end
    end
    qv[i][0] = q
    nextStateTime[i] = simt + h
    if (q-x)>1*Δ
        q=x+1*Δ
    elseif (q-x)<-1*Δ
        q=x-1*Δ
        
    end
    return nothing
end
function Liqss_reComputeNextTime(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64})
    dt = 0.0
    q = qv[i][0]
    x = xv[i][0]
    x1 = xv[i][1]
    if abs(q - x) >= 2 * quantum[i]
        nextStateTime[i] = simt + 1e-12
    else
        if x1 != 0.0
            dt = (q - x) / x1
            if dt > 0.0
                nextStateTime[i] = simt + dt
            elseif dt < 0.0
                if x1 > 0.0
                    nextStateTime[i] = simt + (q - x + 2 * quantum[i]) / x1
                else
                    nextStateTime[i] = simt + (q - x - 2 * quantum[i]) / x1
                end
            end
        else
            nextStateTime[i] = Inf
        end
    end
    if nextStateTime[i] <= simt
        nextStateTime[i] = simt + Inf
    end
end
