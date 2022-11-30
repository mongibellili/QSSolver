#= function updateQ(::Val{1},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    qaux[i][1]=qv[i][0]# index shift....sorry but be careful: taylor 1st elemtn is at 0, a vect 1st elemnt is at 1
    olddx[i][1]=xv[i][1]
   # qv[i][0]=xv[i][0]
    q=xv[i][0]
    u=uv[i][i][1]  # for order 2: u=u+tu*deru
    #dq=0.0
   # x=xv[i][0]
    x1=xv[i][1]
    a=av[i][i]
 #   if a[i][i] !=0.0
        if x1>0
            dx=u+(q+quantum[i])*a
            if dx>=0
                dq=quantum[i]
            else
                dq=(-u/a)-q
            end
        else
            dx=u+(q-quantum[i])*a
            if dx<=0
                dq=-quantum[i]
            else
                dq=(-u/a)-q
            end
        end
    #=     else   #if a==0 u same as derx meaning that in updateQ if derx> we dont have to check if u>0 ....note1
        if x[i][1]>0
            dq=quantum[i]
        else
            dq=-quantum[i]
        end
        end 
    end=#

    if dq>2*quantum[i]
        dq=2*quantum[i]# why not just quantum?
    end
    if dq<-2*quantum[i]
        dq=-2*quantum[i]
    end
    qv[i][0]=q+dq
    return nothing
end =#
function updateQ(::Val{1},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    qaux[i][1]=qv[i][0]# index shift....sorry but be careful: taylor 1st elemtn is at 0, a vect 1st elemnt is at 1
    olddx[i][1]=xv[i][1]
    #q[i][0]=x[i][0]
    q=qv[i][0]
    u=uv[i][i][1]  # for order 2: u=u+tu*deru
    #dq=0.0
    x=xv[i][0]
    dx=xv[i][1]
    a=av[i][i]
    if a !=0.0
        if dx==0.0
           # dx=u+(q)*a
           # if dx==0.0
                dx=1e-26
           # end
        end
    #for order1 finding h is easy but for higher orders iterations are cheaper than finding exact h using a quadratic,cubic...
    #exacte for order1: h=-2Δ/(u+xa-2aΔ) or h=2Δ/(u+xa+2aΔ)
        h = ft-simt
        q = (x + h * u) /(1 - h * a)
        if (abs(q - x) > 2 * quantum[i]) # removing this did nothing...check @btime later
          h = (abs(2 * quantum[i] / dx));
          q= (x + h * u) /(1 - h * a)
        end
        while (abs(q - x) > 2 * quantum[i]) 
          h = h * (quantum[i] / abs(q - x));
          q= (x + h * u) /(1 - h * a)
        end
 
    else
        dx=u
        if dx>0.0
            q=x+quantum[i]# 
        else
            q=x-quantum[i]
        end
    end
    qv[i][0]=q
   # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
    return nothing
end
function updateQ(::Val{2},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent
    tq[i]=simt
    olddx[i][1]=x1#appears only here...updated here and used in updateApprox   
    u1=u1+(simt-tu[i])*u2 # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
    uv[i][i][1]=u1
    tu[i]=simt  
    ddx=x2
    a=av[i][i]
    quan=quantum[i]
    if a!=0.0
        if ddx ==0.0
            #= ddx=a*a*q+a*u1 +u2
            if ddx==0.0 =#
                ddx=1e-26# changing -40 to -6 nothing changed
                println("ddx=0")
            #end
        end
        h = ft-simt
        q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                 (1 - h * a + h * h * a * a / 2)
        if (abs(q - x) >  quan) # removing this did nothing...check @btime later
          h = sqrt(abs(2*quan / ddx)) # sqrt highly recommended...removing it leads to many sim steps..//2* is necessary in 2*quan when using ddx
          q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2)
        end
        while (abs(q - x) >  quan) 
          h = h *quan / abs(q - x)
          q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2)
        end
        q1=(a*q+u1+h*u2)/(1-h*a)  #later investigate 1=h*a
    else
        #ddx=u2
        #= if x2>0.0  #if ddx>0.0 same results
            q=x+quan
        else# elseif x2<0   ......else q=x??
            q=x-quan
        end =#
        if x2!=0.0
            h=sqrt(abs(2*quan/x2))   #sqrt necessary with u2
           # q=x+h*x1+h*h*x2/2
           q=x-h*h*x2/2
            q1=x1+h*x2  #(250 allocations: 16.07 KiB)
        else
            println("x2==0")
            if x1!=0.0
                h=abs(quan/x1)
                q=x+h*x1
                q1=x1
            else
                q=x
                q1=x1
            end
        #=q1=x1+h*x2  #(252 allocations: 16.07 KiB)
        else
            q1=x1 =#
            #println("ddx=0")
        end 
    end
    #= if abs(q-x)>2*quantum[i]#uncommenting this did nothing
        q=x
    end =#
    #olddx[i][2]=ddx  #olddx[i][2] never used again so no need to update it 
    qv[i][0]=q
    qv[i][1]=q1  
   # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
    return nothing
end

#= function updateQ(::Val{2},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    q=xv[i][0];q1=qv[i][1];x=xv[i][0];x1=xv[i][1];x2=xv[i][2]*2;u=uv[i][i][1];u1=uv[i][i][2]
    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent
    tq[i]=simt
    olddx[i][1]=x1#appears only here...updated here and used in updateApprox   
    u=u+(simt-tu[i])*u1 # for order 2: u=u+tu*deru  this is necessary. deleting causes scheduler error
    uv[i][i][1]=u
    tu[i]=simt  
   # olddx[i][2]=2*x2# 
    ddx=x2
    dq=0.0
    a=av[i][i]
    qun=quantum[i]
    if a!=0.0       
        #= if x2>0
            ddx=a*a*(q+qun)+a*u +u1
            if ddx>0.0
                dq=qun# 
            else
                dq=-u1/(a*a)-u/(a)-q
            end
            if abs(dq)>qun
                dq=qun
            end
        elseif x2<0
            ddx=a*a*(q-qun)+a*u +u1
            if ddx<0.0
                dq=-qun# 
            else
                dq=-u1/(a*a)-u/(a)-q
            end
            if abs(dq)>qun
                dq=-qun
            end
        else
            dq=-u1/(a*a)-u/(a)-q
            if abs(dq)>qun
                dq=-qun
            end
        end     =#
        if x2>0
            ddx=a*a*(q+qun)+a*u +u1
            if ddx>=0.0
                dq=qun# 
            else
                dq=-u1/(a*a)-u/(a)-q
            end
            if abs(dq)>qun
                dq=qun
            end
        else
            ddx=a*a*(q-qun)+a*u +u1
            if ddx<=0.0
                dq=-qun# 
            else
                dq=-u1/(a*a)-u/(a)-q
            end
            if abs(dq)>qun
                dq=-qun
            end
        #= else
            dq=-u1/(a*a)-u/(a)-q
            if abs(dq)>qun
                dq=-qun
            end =#
        end    
        #q=q+dq  
        q1=a*q+u  # this is missing +h*x2=+h(a*q1+u1)...as if there is an order drop here
        #q=x+dq # in this approach new q has to be after q1...ie q1 cannot be calculated using new q....infinite loop..program does not terminate
    else       
        #if x2>0.0
        ddx=u1
        if x2>0.0
            dq=-qun
        else
            dq=qun
        end
       # q=q+dq
       # q=x+quantum[i]  #works fine !!!
       # q=x-quantum[i] #errors ...solution escapes up...later test if this behavior is specific to this problem or it is a general thing
       # q=x+2*quantum[i]  #2*Δ errors solution escapes down
       # q=x+quantum[i]   #removing it errors
       if x2!=0.0
            h=sqrt(abs(4*quantum[i]/(x2)))
            q1=u+h*u1
       else
            q1=u
       # println("ddx=0")
       end 
    end
    #= if abs(q-x)>2*quantum[i]#uncommenting this did nothing
        q=x
    end =#
    #q+=dq
    if dq>2*quantum[i]
        dq=2*quantum[i]# why not just quantum?
    end
    if dq<-2*quantum[i]
        dq=-2*quantum[i]
    end
    q+=dq
   # olddx[i][2]=ddx  #olddx[i][2] never used again so no need to update it 
    qv[i][0]=q
    qv[i][1]=q1  
    return nothing
end =#

#= 
function updateQ(::Val{3},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    q=xv[i][0]  # q=x used in liqss1 and liqss3
    q1=qv[i][1]
    q2=qv[i][2]
    #----------------
    x=xv[i][0]
    x1=xv[i][1]
    x2=xv[i][2]*2
    x3=xv[i][3]#*6
    #----------------
    u=uv[i][i][1]
    u1=uv[i][i][2]
    u2=uv[i][i][3]
    #----------------
    qaux[i][1]=q+(simt-tq[i])*q1+(simt-tq[i])*(simt-tq[i])*q2#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1+(simt-tq[i])*q2                     #never used
    qaux[i][3]=q2                                     #never used
    tq[i]=simt
    olddx[i][1]=x1 #olddx[i][1] appears only here...updated here and used in updateApprox
    
    u=u+(simt-tu[i])*u1+(simt-tu[i])*(simt-tu[i])*u2/2  # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
    uv[i][i][1]=u
    u1=u+(simt-tu[i])*u2 #2*u2# if I make a mistake the error is large at begining but it is extremely low afterwards !!!!!!!!
    uv[i][i][2]=u1
    tu[i]=simt
    
   # olddx[i][2]=2*x2#  olddx[i][2] not used elsewhere....so not needed
   # ddx=olddx[i][2]
   dq=0.0
    a=av[i][i]
    if a!=0.0
        if x3>0
            dddx=a*a*a*(q+quantum[i])+a*a*u+a*u1+u2 #*2
            if dddx >=0 
                dq=quantum[i]
            else
                dq=-u2/(a*a*a)-u1/(a*a)-u/a-q  #*2
            end
            if abs(dq)>quantum[i]
                dq=quantum[i]
            end
        else
            dddx=a*a*a*(q-quantum[i])+a*a*u+a*u1+u2  #*2
            if dddx <=0 
                dq=-quantum[i]
            else
                dq=-u2/(a*a*a)-u1/(a*a)-u/a-q  #*2
            end
            if abs(dq)>quantum[i]
                dq=-quantum[i]
            end
        end
        q1=a*q+u
        q2=a*qv[i][1]+u1
        #= if q1*x1 <0.0 && !flag3
            dq=qaux[i][1]-q-abs(quantum[i])*0.1
        else
            dq=qaux[i][1]-q+abs(quantum[i])*0.1
        end =#
    else #a==0
        if x3>0.0
            dq=-quantum[i]
        else
            dq=+quantum[i]
        end
        if x3!=0.0
            h=cbrt(abs(2*quantum[i]/(x3)))
            q1=u+h*u1+h*h*u2/2   #*2
            q2=u1+h*u2  #*2
       else
           #=  q1=x1#u
            q2=x2#u1 =#
            q1=u
            q2=u1
       # println("ddx=0")
       end 



    end
    if dq>2*quantum[i]
        dq=2*quantum[i]# why not just quantum?
    end
    if dq<-2*quantum[i]
        dq=-2*quantum[i]
    end
    q+=dq




    #olddx[i][2]=ddx   #olddx[i][2] never used again so no need to update it 
    qv[i][0]=q
    qv[i][1]=q1 
    qv[i][2]=q2/2   #/2???

    return nothing
end =#

function updateQ(::Val{3},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    q=qv[i][0];q1=qv[i][1];q2=2*qv[i][2];x=xv[i][0];x1=xv[i][1];x2=2*xv[i][2];x3=6*xv[i][3];u1=uv[i][i][1];u2=uv[i][i][2];u3=uv[i][i][3]
    elapsed=simt-tq[i]
    qaux[i][1]=q+elapsed*q1+elapsed*elapsed*q2/2#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1+elapsed*q2   ;qaux[i][3]=q2     #never used
    olddx[i][1]=x1#+(simt-tq[i])*x2 #later change to tx....appears only here...updated here and used in updateApprox   
    tq[i]=simt
    elapsed=simt-tu[i]
    u1=u1+elapsed*u2+elapsed*elapsed*u3/2  
    uv[i][i][1]=u1
    u2=u2+elapsed*u3 #2*u2
    uv[i][i][2]=u2
    tu[i]=simt
    
   # olddx[i][2]=2*x2# 
    dddx=x3
    a=av[i][i]
    quan=quantum[i]
   # println("before q update",abs(q - x) > 2 * quan)
     if a!=0.0
        if dddx ==0.0
           # dddx=a*a*a*(q)+a*a*u+a*u1+u2 #*2
           # if dddx==0.0
                dddx=1e-26# changing -40 to -6 nothing changed
                println("dddx=0")
           # end
        end
        h = ft-simt
        α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
        β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
        q = β/γ
        if (abs(q - x) >  quan) # removing this did nothing...check @btime later
            h = cbrt(abs((6*quan) / dddx));
            α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
            β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
            γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
            q = β/γ
        end
        maxIter=215
        while (abs(q - x) >  quan) && (maxIter>0)
          maxIter-=1
          h = h *quan / abs(q - x);
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
          γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
          q = β/γ
        end
         #= if maxIter < 215
             println("maxiter of mpudate= ",maxIter)
         end =#
        q1=(a*(1-h*a)*q+u1*(1-h*a)-h*h*(a*u2+u3)/2)/(1-h*a+h*h*a*a/2)
        q2=(a*q1+u2+h*u3)/(1-h*a)
    else
        if x3!=0.0
            h=cbrt(abs(6*quan/x3))
            q=x-h*h*h*x3/6
            q1=x1+x3*h*h/2   #*2
            q2=x2+h*x3
           #=  q=x+h*h*h*u3/6
            q1=u1-h*h*u3/2   
            q2=u2+h*u3   =#
            #= q2=u2+h*u3 
            q1=u1+h*u2+h*h*u3/2-h*q2
            q=x+h*u1+h*h*u2/2+h*h*h*u3/6-h*q1-h*h*q2/2 =#
          #  q=x+x1*h+x2*h*h/2+h*h*h*x3/6        
       else
           #=  q1=x1#u
            q2=x2#u1 =#
            println("x3=0")
            if x2!=0
            
            h=sqrt(abs(2*quan/x2))
            q=x-h*h*x2/2
            q1=x1+h*x2
            
            else
                if x1!=0
                    h=abs(quan/x1)
                    q=x+h*x1
                else
                    q=x
                    
                end
                q1=x1
                
            end
            q2=x2
            # println("ddx=0")
        end 
        #dddx=u2
        #= if x3>0.0
            q=x-quan
        else
            q=x+quan
        end
       # q=x+quantum[i]  #works fine !!!
       # q=x-quantum[i] #errors ...solution escapes up...later test if this behavior is specific to this problem or it is a general thing
       # q=x+2*quantum[i]  #2*Δ errors solution escapes down
       # q=x+quantum[i]   #removing it errors
        if x3!=0.0
            h=cbrt(abs(quan/x3))
            q1=x1+h*x2+h*h*x3/2   #*2
            q2=x2+h*x3  #*2
       else
           #=  q1=x1#u
            q2=x2#u1 =#
            q1=x1
            q2=x2
       # println("ddx=0")
       end  =#
    end
    #= if abs(q-x)>2*quantum[i]#uncommenting this did nothing
        q=x
    end =#
    #olddx[i][2]=ddx  #olddx[i][2] never used again so no need to update it 
    qv[i][0]=q
    qv[i][1]=q1 
    qv[i][2]=q2/2  
    return nothing
end


function Liqss_ComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
    
   
    coef=@SVector [q - x , q1-x1,-x2]#
    nextTime[i] = currentTime + minPosRoot(coef, Val(2))
    
    #nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end
function Liqss_ComputeNextTime(::Val{3}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
    q2=qv[i][2]
    x3=xv[i][3]
   
    coef=@SVector [q - x , q1-x1,q2-x2,-x3]#
    nextTime[i] = currentTime + minPosRoot(coef, Val(3))
    
    #nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end





##########################################################################################################################################################
function Liqss_reComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    dt=0.0
    q=qv[i][0]
    x=xv[i][0]
    if xv[i][1] !=0.0
        dt=(q-x)/xv[i][1]
        if dt>0.0
            nextTime[i]=currentTime+dt# later guard against very small dt
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

#this shortcut did not work (infinite loop) for system15 ft=160.

function Liqss_reComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
    #if xv[i][1] !=0.0
    if  x2!=0.0 && a[i][i] != 0.0  # a!=0 cuz when a==0 x1=0+u  and then u updated below as u=x1-0 meaning that x1==u=oldx1....no change ie nexttime should be Inf???
        nextTime[i]=currentTime+sqrt(abs((q-x)/(x2)))  #64.829 μs (580 allocations: 46.54 KiB)
       #nextTime[i]=currentTime+sqrt(abs((quantum[i])/(x2)))  #64.503 μs (584 allocations: 46.66 KiB)
       #nextTime[i]=currentTime+(abs((x1)/(x2)))#nextTime[i]=currentTime+sqrt(abs((x1)/(x2))) # 54.390 μs (584 allocations: 46.66 KiB)
      # nextTime[i]=currentTime+(abs((q1-x1)/(x2)))# 57.187 μs (586 allocations: 46.73 KiB)
       #coef=@SVector [q - x , q1-x1,-x2]#135.287 μs (580 allocations: 46.54 KiB)
       #nextTime[i] = currentTime + minPosRoot(coef, Val(2))
    else
        nextTime[i]=Inf
    end
    #removing the previous 6 lines and adding time3 -->78.725 μs (656 allocations: 48.91 KiB)
    coef=@SVector [q - x + 2*quantum[i], q1-x1,-x2]#
    time1 = currentTime + minPosRoot(coef, Val(2))
    coef=setindex(coef,q - x - 2*quantum[i],1)
    time2 = currentTime + minPosRoot(coef, Val(2))
    time1 = time1 < time2 ? time1 : time2    
   #=coef=setindex(coef,q - x,1)
    time3 = currentTime + minPosRoot(coef, Val(2))
    time1 = time1 < time3 ? time1 : time3  =#  
    nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
   #=  else
        nextTime[i]=Inf
    end =#
end
#= function Liqss_reComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
 
    coef=@SVector [q - x + 2*quantum[i], q1-x1,-x2]#
    time1 = currentTime + minPosRoot(coef, Val(2))
    coef=setindex(coef,q - x - 2*quantum[i],1)
    time2 = currentTime + minPosRoot(coef, Val(2))
    time1 = time1 < time2 ? time1 : time2    
    coef=setindex(coef,q - x,1)
    time3 = currentTime + minPosRoot(coef, Val(2))
    nextTime[i] = time1 < time3 ? time1 : time3    
    #nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]


    if nextTime[i] < currentTime  
        println("q neither heading towards x nor getting away by 2quantums")
        nextTime[i]=currentTime+1e-6
    end



    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end =#


#= function Liqss_reComputeNextTime(::Val{3}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
    q2=qv[i][2]
    x3=xv[i][3]
   
    coef=@SVector [q - x + 2*quantum[i], q1-x1,q2-x2,-x3]#
    time1 = currentTime + minPosRoot(coef, Val(3))
    coef=setindex(coef,q - x - 2*quantum[i],1)
    time2 = currentTime + minPosRoot(coef, Val(3))
    time1 = time1 < time2 ? time1 : time2   
    coef=setindex(coef,q - x,1)
    time3 = currentTime + minPosRoot(coef, Val(3))
    nextTime[i] = time1 < time3 ? time1 : time3    
    #nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    if nextTime[i] < currentTime  
        nextTime[i]=currentTime+1e-6
    end
    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end =#

function Liqss_reComputeNextTime(::Val{3}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
    q2=qv[i][2]
    x3=xv[i][3]
    if  x3!=0.0 && a[i][i] != 0.0
        #nextTime[i]=currentTime+abs((q2-x2)/(x3))#4.247 ms (934 allocations: 64.02 KiB)
      #  nextTime[i]=currentTime+abs((x2)/(x3))#619.586 μs (934 allocations: 64.02 KiB)#shifts down
       # nextTime[i]=currentTime+(abs((q1-x1)/(x3)))#758.059 ms (934 allocations: 64.02 KiB)
      # nextTime[i]=currentTime+sqrt(abs((q1-x1)/(x3)))#2.997 ms (934 allocations: 64.02 KiB)#correct
       #nextTime[i]=currentTime+sqrt(abs((x1)/(x3)))#654.351 μs (934 allocations: 64.02 KiB)#shifts up a little
       nextTime[i]=currentTime+cbrt(abs((q-x)/(x3)))#634.635 μs (934 allocations: 64.02 KiB)#*******************correct*******************
     # nextTime[i]=currentTime+cbrt(abs((quantum[i])/(x3)))#632.743 μs (934 allocations: 64.02 KiB)#shifts up
    else
        nextTime[i]=Inf
    end
    coef=@SVector [q - x + 2*quantum[i], q1-x1,q2-x2,-x3]#
    time1 = currentTime + minPosRoot(coef, Val(3))
    coef=setindex(coef,q - x - 2*quantum[i],1)
    time2 = currentTime + minPosRoot(coef, Val(3))
    time1 = time1 < time2 ? time1 : time2    
    nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end


#######################################################################################################################################################
function updateLinearApprox(::Val{1},i::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=q[i][0]-qaux[i][1]
    if diffQ != 0.0
        a[i][i]=(x[i][1]-olddx[i][1])/diffQ
    else
        a[i][i]=0.0
    end
    u[i][i][1]=x[i][1]-a[i][i]*q[i][0]   #if a==0 u same as derx meaning that in updateQ if derx> we dont have to check if u>0 ....note1
    return nothing
end
function updateLinearApprox(::Val{2},i::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=q[i][0]-qaux[i][1]
   # println("aii before updateoher= ",a[i][i])
    if diffQ != 0.0
        a[i][i]=(x[i][1]-olddx[i][1])/diffQ
       #=  if a[i][i]>0.0   # this is new for liqss2....deleting this causes the schedule static error unless changing also inside recompute a<0 to a!=0
            a[i][i]=0.0
        end =#
        #= if a[i][i]==0.0   # this is new for liqss2....deleting this causes the schedule static error unless changing also inside recompute a<0 to a!=0
            println("x[i][1]==olddx[i][1] at simt= ",simt)
        end =#
    else
        a[i][i]=0.0
       # println("diffQ==0 at simt= ",simt)
        #x[i][1]=0.0
        
    end
   # println("aii before updateoher= ",a[i][i])
    u[i][i][1]=x[i][1]-a[i][i]*q[i][0]
    u[i][i][2]=2*x[i][2]-a[i][i]*q[i][1]
    #tu[i]=simt  # comment did nothing but it makes sense to keep it because more accurate since u is changed
    if debug
        println("u$i$i = ",u[i][i][1])
        println("du$i$i = ",u[i][i][2])
        println("tu[$i] = ",tu[i])
    end
    return nothing
end
function updateLinearApprox(::Val{3},i::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=q[i][0]-qaux[i][1]
    if diffQ != 0.0  #if (fabs(diffQ) > lqu[var] * 1e-6)
        a[i][i]=(x[i][1]-olddx[i][1])/diffQ
       #=  if a[i][i]>0.0   # this is new for liqss2....deleting this causes the schedule static error unless changing also inside recompute a<0 to a!=0
            a[i][i]=0.0
        end =#
        #= if a[i][i]==0.0   # this is new for liqss2....deleting this causes the schedule static error unless changing also inside recompute a<0 to a!=0
            println("x[i][1]==olddx[i][1] at simt= ",simt)
        end =#
    else
        a[i][i]=0.0
       # println("diffQ==0 at simt= ",simt)
       # x[i][1]=0.0
    end
    u[i][i][1]=x[i][1]-a[i][i]*q[i][0]    
    u[i][i][2]=2*x[i][2]-a[i][i]*q[i][1]  #2*x[i][2]-a[i][i]*q[i][1]  #if a==0 deru same as derderx meaning that in updateQ if derderx> we dont have to check if deru>0 ....note1
    u[i][i][3]=6*x[i][3]-a[i][i]*2*q[i][2]  #3*x[i][3]-a[i][i]*q[i][2]
    #tu[i]=simt  # uncomment did nothing
   
    return nothing
end
