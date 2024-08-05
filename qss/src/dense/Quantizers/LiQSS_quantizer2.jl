#analytic
function updateQ(::Val{2},::Val{1}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})
    error("a quick check for future derivative does not exist in order 2 updateQ")
    return nothing
end
function updateQ(::Val{2},::Val{2}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})
    cacheA[1] = 0.0
    exactA(qv, d, cacheA, i, i, simt)
    a = cacheA[1]
    # exactA(xv,cacheA,i,i);a=cacheA[1]
    q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; #u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent

    u1=x1-a*qaux[i][1]
    u2=x2-a*q1
    dxaux[i][1]=x1
    dxaux[i][2]=x2

    ddx=x2
    quan=quantum[i]
    h=0.0
    if a!=0.0
        if ddx ==0.0
            ddx=a*a*q+a*u1 +u2
            if ddx==0.0 
                ddx=1e-40# changing -40 to -6 nothing changed
                # println("ddx=0")
            end
        end
        h = ft-simt; q=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
                
        
        if (abs(q - x) > 2.0* quan) # removing this did nothing...check @btime later
        h = sqrt(abs(2*quan / ddx)) # sqrt highly recommended...removing it leads to many sim steps..//2* is necessary in 2*quan when using ddx
        #= q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                    (1 - h * a + h * h * a * a / 2) =#
                #   qtemp=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
                    q=x+(-h*h*(a*a*x+a*u1+u2)/2)/((1 - h * a + h * h * a * a / 2))
                #=  if abs(q-qtemp)>1e-16
                    @show qtemp,q
                    end =#
                
        end
        maxIter=10000
        # tempH=h
        while (abs(q - x) >2.0*  quan) && (maxIter>0) && (h>0)
            
        h = h *sqrt(quan / abs(q - x))
        #  h = h *0.99*(2.0*quan / abs(q - x))
        #= q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                    (1 - h * a + h * h * a * a / 2) =#
                #  qtemp=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
                    q=x+(-h*h*(a*a*x+a*u1+u2)/2)/((1 - h * a + h * h * a * a / 2))
                    #= if abs(q-qtemp)>1e-16
                    @show q,qtemp
                    end =#
                
        maxIter-=1
        end
        if maxIter==0
            println("updateQ maxIterReached")
        end
    
            

        

        α1=1-h*a
        if abs(α1)==0.0
            α1=1e-30*sign(α1)
        end
        q1=(a*q+u1+h*u2)/α1  #later investigate 1=h*a


    else
        println("a==0")
        if x2!=0.0  
            h=sqrt(abs(2*quan/x2))   #sqrt necessary with u2
            q=x-h*h*x2/2
            q1=x1+h*x2
        else
            # println("x2==0")
            if x1!=0.0
                h=abs(quan/x1)
                q=x+h*x1
                q1=0#x#/2
            else
                h=Inf
                q=x
                q1=x1
            end
        end 

    end

    qv[i][0]=q
    qv[i][1]=q1  
    nextStateTime[i]=simt+h

 
    return nothing
end
function updateQ(::Val{2},::Val{3}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})
    cacheA[1] = 0.0
    exactA(qv, d, cacheA, i, i, simt)
    a = cacheA[1]
    # exactA(xv,cacheA,i,i);a=cacheA[1]
    q = qv[i][0]
    q1 = qv[i][1]
    x = xv[i][0]
    x1 = xv[i][1]
    x2 = xv[i][2] * 2 #u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1] = q + (simt - tq[i]) * q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2] = q1                     #appears only here...updated here and used in updateQevent
    u1 = x1 - a * qaux[i][1]
    u2 = x2 - a * q1
    dxaux[i][1] = x1
    dxaux[i][2] = x2
    ddx = x2
    quan = quantum[i]
    h = 0.0
    if a != 0.0
        if a * a * x + a * u1 + u2 <= 0.0
            if -(a * a * x + a * u1 + u2) / (a * a) < quan # asymptote<delta...no sol ...no need to check
                h = Inf
                q = -(a * u1 + u2) / (a * a)#q=x+asymp
            else
                q = x + quan
                coefi = NTuple{3,Float64}(((a * a * (x + quan) + a * u1 + u2) / 2, -a * quan, quan))
                h = minPosRootv1(coefi) #needed for dq
                #math was used to show h cannot be 1/a
            end
        elseif a * a * x + a * u1 + u2 >= 0.0
            if -(a * a * x + a * u1 + u2) / (a * a) > -quan # asymptote>-delta...no sol ...no need to check
                h = Inf
                q = -(a * u1 + u2) / (a * a)
            else
                coefi = NTuple{3,Float64}(((a * a * (x - quan) + a * u1 + u2) / 2, a * quan, -quan))
                h = minPosRootv1(coefi)
                q = x - quan
            end
        else#a*a*x+a*u1+u2==0 -->f=0....q=x+f(h)=x+h*g(h) -->g(h)==0...dx==0 -->aq+u==0
            q = -u1 / a
            h = Inf
        end
        if h != Inf
            if h != 1 / a
                q1 = (a * q + u1 + h * u2) / (1 - h * a)
            else#h=1/a
                error("report bug: updateQ: h cannot=1/a; h=$h , 1/a=$(1/h)")
            end
        else #h==inf make ddx==0 dq=-u2/a
            q1 = -u2 / a
        end
    else
        #println("a==0")
        if x2 != 0.0
            h = sqrt(abs(2 * quan / x2))   #sqrt necessary with u2
            q = x - h * h * x2 / 2
            q1 = x1 + h * x2
        else
            # println("x2==0")
            if x1 != 0.0
                #quantum[i]=1quan
                h = abs(1 * quan / x1)   # *10 just to widen the step otherwise it would behave like 1st order
                q = x + h * x1
                q1 = 0.0
            else
                h = Inf
                q = x
                q1 = 0.0
            end
        end
    end
    qv[i][0] = q
    qv[i][1] = q1
    nextStateTime[i] = simt + h
    return h
end
function Liqss_reComputeNextTime(::Val{2}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}) #= ,a::Vector{Vector{Float64}} =#
    q = qv[i][0]
    x = xv[i][0]
    q1 = qv[i][1]
    x1 = xv[i][1]
    x2 = xv[i][2]
    quani = quantum[i]
    β = 0
    if abs(q - x) >= 2 * quani # this happened when var i and j s turns are now...var i depends on j, j is asked here for next time...or if you want to increase quant*10 later it can be put back to normal and q & x are spread out by 10quan
        nextStateTime[i] = simt + 1e-12
    else
        coef = @SVector [q - x, q1 - x1, -x2]#
        nextStateTime[i] = simt + minPosRoot(coef, Val(2))
        if q - x > 0.0#1e-9
            coef = setindex(coef, q - x - 2 * quantum[i], 1)
            timetemp = simt + minPosRoot(coef, Val(2))
            if timetemp < nextStateTime[i]
                nextStateTime[i] = timetemp
            end
        elseif q - x < 0.0#-1e-9
            coef = setindex(coef, q - x + 2 * quantum[i], 1)
            timetemp = simt + minPosRoot(coef, Val(2))
            if timetemp < nextStateTime[i]
                nextStateTime[i] = timetemp
            end
        end
        if nextStateTime[i] <= simt # this is coming from the fact that a variable can reach 2quan distance when it is not its turn, then computation above gives next=simt+(p-p)/dx...p-p should be zero but it can be very small negative
            nextStateTime[i] = simt + Inf#1e-14
            # @show simt,nextStateTime[i],i,x,q,quantum[i],xv[i][1]
        end
    end
end
