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
    if a != 0.0
        qplus=x+sign(x1)*Δ
        dxplus=a*(qplus)+u
    
        if x1*dxplus>0.0
            q=qplus
            h=Δ/abs(dxplus)
        else
            q=-u/a
            h=Inf
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
    if (q-x)>1*Δ # this is for q=-u/a
        q=x+1*Δ
    elseif (q-x)<-1*Δ
        q=x-1*Δ
    end
    qv[i][0] = q
    nextStateTime[i] = simt + h
    return nothing
end






function updateQ(::Val{1},opt::Options{2,DU}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64},a::Float64, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})where{DU}
    
    q = qv[i][0]
    x = xv[i][0]
    x1 = xv[i][1]
    qaux[i][1] = q
    u = x1 - a * q
    dxaux[i][1] = x1
    h = 0.0
    Δ = quantum[i]
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
    if (q-x)>1*Δ # this is for q=-u/a
        q=x+1*Δ
    elseif (q-x)<-1*Δ
        q=x-1*Δ
    end
    qv[i][0] = q
    nextStateTime[i] = simt + h
    return nothing
end


function updateQ(::Val{1},opt::Options{3,DU}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64},a::Float64, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})where{DU}
    q = qv[i][0]
    x = xv[i][0]
    x1 = xv[i][1]
    qaux[i][1] = q
    u = x1 - a * q
    dxaux[i][1] = x1
    dx=x1
    h=0.0
    Δ=quantum[i]
   if a !=0.0
       h = ft-simt
       if (1 - h * a)!=0 q = (x + h * u) /(1 - h * a) end
      if (abs(q - x) > 1* Δ) 
         h = (abs( Δ / x1));
         if (1 - h * a)!=0 q= (x + h * u) /(1 - h * a) end
       end
       if (abs(q - x) > 1* Δ) 
         h = h *sqrt(Δ / abs(q - x))
         if (1 - h * a)!=0 q= (x + h * u) /(1 - h * a)  end
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
    if (q-x)>1*Δ# this is for "if (abs(q - x) > 1* Δ) "
        q=x+1*Δ
      elseif (q-x)<-1*Δ
        q=x-1*Δ
      end
    qv[i][0] = q
    nextStateTime[i] = simt + h
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
