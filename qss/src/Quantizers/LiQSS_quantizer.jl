function updateQ(::Val{1},i::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{O,Float64}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    qaux[i][1]=q[i][0]# index shift....sorry but be careful: taylor 1st elemtn is at 0, a vect 1st elemnt is at 1
    olddx[i][1]=x[i][1]
    q[i][0]=x[i][0]
    q_=q[i][0]
    u_=u[i][1]  # for order 2: u=u+tu*deru
    dq=0.0
    if x[i][1]>0
        dx=u_+(q_+quantum[i])*a[i][i]
        if dx>=0
            dq=quantum[i]
        else
            dq=(-u_/a[i][i])-q_
        end
    else
        dx=u_+(q_-quantum[i])*a[i][i]
        if dx<=0
            dq=-quantum[i]
        else
            dq=(-u_/a[i][i])-q_
        end
    end
    if dq>2*quantum[i]
        dq=2*quantum[i]# why not just quantum?
    end
    if dq<-2*quantum[i]
        dq=-2*quantum[i]
    end
    q[i][0]=q_+dq
    return nothing
end
function updateQ(::Val{2},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{O,Float64}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
   
    q=qv[i][0]
    q1=qv[i][1]
    x=xv[i][0]
    x1=xv[i][1]
    x2=xv[i][2]
    u=uv[i][1]
    u1=uv[i][2]

    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent
    tq[i]=simt
    olddx[i][1]=x1#appears only here...updated here and used in updateApprox
    
    u=u+(simt-tu[i])*u1 # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
    uv[i][1]=u
    tu[i]=simt
    
    olddx[i][2]=2*xv[i][2]# 
    ddx=olddx[i][2]
    a=av[i][i]
     if a!=0.0
        if ddx ==0.0
            ddx=a*a*q+a*u +u1
            if ddx==0.0
                ddx=1e-26# changing -40 to -6 nothing changed
            end
        end
        h=ft-simt
        #q=((x+h*u+h*h*u1/2)*(1-h*a)+(h*h*a/2-h)*(u+h*u1))/(1-h*a+h*h*a*a/2)
        #q=((x+h*u+h*h*u1/2)*(1-h*a)+(h*h*a/2)*(u+h*u1))/(1-2*h*a+h*h*a*a/2)
        q=(x+h*u+h*h*(a*u+u1)/2)/(1-h*a-h*h*a*a/2)
        if abs(q-x)>2*quantum[i]
            h=sqrt(abs(2*quantum[i]/ddx))
            #q=((x+h*u+h*h*u1/2)*(1-h*a)+(h*h*a/2-h)*(u+h*u1))/(1-h*a+h*h*a*a/2)
            #q=((x+h*u+h*h*u1/2)*(1-h*a)+(h*h*a/2)*(u+h*u1))/(1-2*h*a+h*h*a*a/2)
            q=(x+h*u+h*h*(a*u+u1)/2)/(1-h*a-h*h*a*a/2)
        end
#=         @timeit "updatQ-while" while abs(q-x)>2*quantum[i]
            h=h*sqrt(abs(2*quantum[i]/(q-x)))
            #q=((x+h*u+h*h*u1/2)*(1-h*a)+(h*h*a/2-h)*(u+h*u1))/(1-h*a+h*h*a*a/2)
            #q=((x+h*u+h*h*u1/2)*(1-h*a)+(h*h*a/2)*(u+h*u1))/(1-2*h*a+h*h*a*a/2)
            q=(x+h*u+h*h*(a*u+u1)/2)/(1-h*a-h*h*a*a/2)
        end =#
        coef=@SVector [ -2* quantum[i], x*a+2*quantum[i]+u,(a*u+u1+x*a*a+2*quantum[i]*a*a)/2]#*2
        h =  minPosRoot(coef, Val(2))
        if h==Inf
             coef=@SVector [ 2* quantum[i], x*a-2*quantum[i]+u,(a*u+u1+x*a*a-2*quantum[i]*a*a)/2]#*2
             h =  minPosRoot(coef, Val(2))
        end
        if h==Inf  #probably not needed
            while abs(q-x)>2*quantum[i]
                h=h*sqrt(abs(2*quantum[i]/(q-x)))
                #q=((x+h*u+h*h*u1/2)*(1-h*a)+(h*h*a/2-h)*(u+h*u1))/(1-h*a+h*h*a*a/2)
                #q=((x+h*u+h*h*u1/2)*(1-h*a)+(h*h*a/2)*(u+h*u1))/(1-2*h*a+h*h*a*a/2)
                q=(x+h*u+h*h*(a*u+u1)/2)/(1-h*a-h*h*a*a/2)
            end
        end

        q=(x+h*u+h*h*(a*u+u1)/2)/(1-h*a-h*h*a*a/2)
        q1=(a*q+u+h*u1)/(1-h*a)
 #=        if x2*ddx<0.0  #not necessary :)#uncommenting did nothing
            q1=-u1/a
            q=(q1-u)/a
            if abs(q-x)>2*quantum[i]
                if q>x
                    q=x+quantum[i]
                else
                    q=x-quantum[i]
                end
            end
        end =#
         else
      #  println("a=0")
       # println("simt= ",simt)
        ddx=u1
        if ddx>0.0
            q=x-quantum[i]
        else
            q=x+quantum[i]
        end
        if ddx!=0.0
            h=sqrt(abs(2*quantum[i]/ddx))
            q1=u+h*u1
        else
            q1=u
           # println("ddx=0")
        end
    end
    if abs(q-x)>2*quantum[i]#uncommenting this did nothing
        q=x
    end
    olddx[i][2]=ddx
     qv[i][0]=q
    qv[i][1]=q1  
    return nothing
end


##########################################################################################################################################################
function Liqss_reComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    dt=0.0
    q=qv[i][0]
    x=xv[i][0]
    if xv[i][1] !=0.0
        dt=(q-x)/xv[i][1]
        if dt>0.0
            nextTime[i]=currentTime+dt
        else
            if xv[i][1]>0.0  
                nextTime[i]=currentTime+(q-x+2*quantum[i])/xv[i][1]
            else
                nextTime[i]=currentTime+(q-x-2*quantum[i])/xv[i][1]
            end
            if nextTime[i] < currentTime  
                nextTime[i]=currentTime+1e-6
            end
        end
    else
        nextTime[i]=Inf
    end
end

function Liqss_reComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
    if  x2!=0.0 && a[i][i] != 0.0
        nextTime[i]=currentTime+abs((q1-x1)/x2)
    else
        nextTime[i]=Inf
    end
    coef=@SVector [q - x + 2*quantum[i], q1-x1,-x2]#*2 or not nothing changes
    time1 = currentTime + minPosRoot(coef, Val(2))
    coef=setindex(coef,q - x - 2*quantum[i],1)
    time2 = currentTime + minPosRoot(coef, Val(2))
    time1 = time1 < time2 ? time1 : time2    
    nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end
#######################################################################################################################################################
function updateLinearApprox(::Val{1},i::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{O,Float64}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=q[i][0]-qaux[i][1]
    if diffQ != 0.0
        a[i][i]=(x[i][1]-olddx[i][1])/diffQ
    else
        a[i][i]=0.0
    end
    u[i][1]=x[i][1]-a[i][i]*q[i][0]
    return nothing
end

function updateLinearApprox(::Val{2},i::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{O,Float64}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=q[i][0]-qaux[i][1]
    if diffQ != 0.0
        a[i][i]=(x[i][1]-olddx[i][1])/diffQ
       #=  if a[i][i]>0.0   # this is new for liqss2....deleting this causes the schedule static error unless changing also inside recompute a<0 to a!=0
            a[i][i]=0.0
        end =#
    else
        a[i][i]=0.0
    end
    u[i][1]=x[i][1]-a[i][i]*q[i][0]
    u[i][2]=x[i][2]-a[i][i]*q[i][1]
   # tu[i]=simt  # uncomment did nothing
    return nothing
end