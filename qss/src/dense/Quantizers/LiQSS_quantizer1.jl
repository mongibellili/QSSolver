

function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
        cacheA[1]=0.0;exacteA(qv,cacheA,i,i);a=cacheA[1]
         q=qv[i][0];x=xv[i][0];x1=xv[i][1];
         qaux[i][1]=q
         u=x1-a*q
         dxaux[i][1]=x1
         h=0.0
         Δ=quantum[i]
    
    u=x1-a*q
    #uv[i][i][1]=u
   
    dxaux[i][1]=x1
   qplus=x+sign(x1)*quantum[i]
   dxplus=a*(qplus)+u
   if x1*dxplus>0.0
         qv[i][0]=qplus
    else
        qv[i][0]=-u/a
    end
    return nothing
end
 

 
function updateQ(::Val{2},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    cacheA[1]=0.0;exacteA(qv,cacheA,i,i);a=cacheA[1]
     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
     u=x1-a*q
     dxaux[i][1]=x1
     dx=x1
     h=0.0
     Δ=quantum[i]
     
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
       if (abs(q - x) > 1* quantum[i]) 
          h = (abs( quantum[i] / dx));
          q= (x + h * u) /(1 - h * a)
        end
        while (abs(q - x) > 1* quantum[i]) 
          h = h * 0.98*(quantum[i] / abs(q - x));
          q= (x + h * u) /(1 - h * a)
        end 
       
    else
        dx=u
        if dx>0.0
            q=x+quantum[i]# 
        else
            q=x-quantum[i]
        end
        if dx!=0
        h=(abs(quantum[i]/dx))
        else
            h=Inf
        end
    end
    qv[i][0]=q
   # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
   nextStateTime[i]=simt+h
    return nothing
end

function updateQ(::Val{3},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    cacheA[1]=0.0;exacteA(qv,cacheA,i,i);a=cacheA[1]
     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
     u=x1-a*q
     dxaux[i][1]=x1
     h=0.0
     Δ=quantum[i]
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
    qv[i][0]=q
   # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
   nextStateTime[i]=simt+h
    return nothing
end



  





















function Liqss_ComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64})
    q=qv[i][0];x=xv[i][0];x1=xv[i][1]
    if  x1!=0.0  
        nextStateTime[i]=currentTime+(abs((q-x)/(x1)))  
    else
        nextStateTime[i]=Inf
    end
end

function Liqss_reComputeNextTime(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,a::Vector{Vector{Float64}} =#)
    dt=0.0; q=qv[i][0];x=xv[i][0];x1=xv[i][1]
    if abs(q-x) >= 2*quantum[i] # this happened when var i and j s turns are now...var i depends on j, j is asked here for next time...or if you want to increase quant*10 later it can be put back to normal and q & x are spread out by 10quan
        nextStateTime[i] = simt+1e-14
    else
                if x1 !=0.0 #&& abs(q-x)>quantum[i]/10
                    dt=(q-x)/x1
                    if dt>0.0
                        nextStateTime[i]=simt+dt# later guard against very small dt
                    elseif dt<0.0
                        if x1>0.0  
                            nextStateTime[i]=simt+(q-x+2*quantum[i])/x1
                        else
                            nextStateTime[i]=simt+(q-x-2*quantum[i])/x1
                        end
                    end
                else
                    nextStateTime[i]=Inf
                end
    end
  #=   if nextStateTime[i]<simt # unnecessary never reached
        nextStateTime[i]=simt+Inf#1e-12
    end =#
  # this is coming from the fact that a variable can reach 2quan distance when it is not its turn, then computation above gives next=simt+(p-p)/dx...p-p should be zero but it can be very small negative
end



