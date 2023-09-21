 #iters
#=  function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
   # a=av[i][i]
     exacteA(qv,cacheA,i,i)
    
     a=cacheA[1]
 

    #=    if i==1
        if abs(a+1.1)>1e-3
            @show i,a
        end
    else
        if abs(a+20.0)>1e-3
            @show i,a
        end
    end =#

     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
    #=  olddx[i][1]=x1 =#
     u=x1-a*q
     #uv[i][i][1]=u
     dx=x1
     dxaux[i][1]=x1
     h=0.0
     Δ=quantum[i]
     debugH=4.3;debugL=14.1
    if a !=0.0
        if dx==0.0
            dx=u+(q)*a
            if dx==0.0
                dx=1e-26
            end
        end
    #for order1 finding h is easy but for higher orders iterations are cheaper than finding exact h using a quadratic,cubic...
    #exacte for order1: h=-2Δ/(u+xa-2aΔ) or h=2Δ/(u+xa+2aΔ)
        h = ft-simt
        q = (x + h * u) /(1 - h * a)
        if (abs(q - x) >  1*quantum[i]) # removing this did nothing...check @btime later
          h = (abs( quantum[i] / dx));
          q= (x + h * u) /(1 - h * a)
        end
        while (abs(q - x) >  1*quantum[i]) 
          h = h * 0.99*(quantum[i] / abs(q - x));
          q= (x + h * u) /(1 - h * a)
        end
       
       #=  α=-(a*x+u)/a
        h1denom=a*(x+Δ)+u;h2denom=a*(x-Δ)+u
        if h1denom==0.0 h1denom=1e-26 end
        if h2denom==0.0 h2denom=1e-26 end
        h1=Δ/h1denom
        h2=-Δ/h2denom
        if a<0
           if α>Δ
               h=h1;q=x+Δ
           elseif α<-Δ
               h=h2;q=x-Δ
           else
               h=max(ft-simt,-1/a)
               q=(x+h*u)/(1-h*a)#1-h*a non neg because h > -1/a > 1/a
           end
        else #a>0
            if α>Δ
                h=h2;q=x-Δ
            elseif α<-Δ
                h=h1;q=x+Δ
            else
                if a*x+u>0
                    h=max(ft-simt,-(a*x+u+a*Δ)/(a*(a*x+u-a*Δ)))# midpoint between asymptote and delta
                else
                    h=max(ft-simt,-(a*x+u-a*Δ)/(a*(a*x+u+a*Δ)))
                end
                q=(x+h*u)/(1-h*a)#1-h*a non neg because values of h from graph show that h >1/a
            end
        end =#
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
    return h
end  =#
 
 
 #analytic favor q-x
   function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
     exacteA(qv,cacheA,i,i);a=cacheA[1]
     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
     u=x1-a*q
     dxaux[i][1]=x1
     h=0.0
     Δ=quantum[i]
     if a !=0.0
         α=-(a*x+u)/a
         h1denom=a*(x+Δ)+u;h2denom=a*(x-Δ)+u
         if h1denom==0.0 h1denom=1e-26 end
         if h2denom==0.0 h2denom=1e-26 end
         h1=Δ/h1denom
         h2=-Δ/h2denom
         if a<0
            if α>Δ
                h=h1;q=x+Δ
            elseif α<-Δ
                h=h2;q=x-Δ
            else
                h=Inf
                q=-u/a
            end
        else #a>0
            if α>Δ
                h=h2;q=x-Δ
            elseif α<-Δ
                h=h1;q=x+Δ
            else
                if a*x+u>0
                    q=x-Δ   #
                    h=h2
                else
                    q=x+Δ 
                    h=h1
                end
            end
        end
     else #a=0
         if x1>0.0
             q=x+Δ # 
         else
             q=x-Δ 
         end
         if x1!=0
            h=(abs(Δ /x1))
         else
             h=Inf
         end
     end
    qv[i][0]=q
    nextStateTime[i]=simt+h
    return h
end   
 
 #analytic favor q-x but care about h large
#= function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    
    # a=av[i][i]
     exacteA(qv,cacheA,i,i)
    
     a=cacheA[1]
 

  #=    if i==1
        if abs(a+1.1)>1e-3
            @show i,a
        end
    else
        if abs(a+20.0)>1e-3
            @show i,a
        end
    end =#

     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
    #=  olddx[i][1]=x1 =#
     u=x1-a*q
     #uv[i][i][1]=u
     dx=x1
     dxaux[i][1]=x1
     h=0.0
     Δ=quantum[i]
     if a !=0.0
         if dx==0.0
             dx=u+(q)*a
             if dx==0.0
                 dx=1e-26
             end
         end
         α=-(a*x+u)/a
         h1denom=a*(x+Δ)+u;h2denom=a*(x-Δ)+u
         if h1denom==0.0 h1denom=1e-26 end
         if h2denom==0.0 h2denom=1e-26 end
         h1=Δ/h1denom
         h2=-Δ/h2denom
         if a<0
            if α>Δ
                h=h1;q=x+Δ
            elseif α<-Δ
                h=h2;q=x-Δ
            else
               #=  h=Inf
                if a*x+u>0
                    q=x+α
                else
                    q=x+α
                end =#
                
                
                   #=  h=-1/a
                    q=(x+h*u)/(1-h*a) =#
                
                   
                h=max(ft-simt,-1/a)
                q=(x+h*u)/(1-h*a)#1-h*a non neg because a<0
                
            end
        else #a>0
            if α>Δ
                h=h2;q=x-Δ
            elseif α<-Δ
                h=h1;q=x+Δ
            else
               #=  h=Inf
                if a*x+u>0
                    q=x+α
                else
                    q=x+α
                end =#
                
                #= if a*x+u>0
                    h=h2;q=x-Δ
                else
                    h=h1;q=x+Δ
                end =#
               #=  h=1/a+ft-simt # i have if simt>=ft break in intgrator
                q=(x+h*u)/(1-h*a) =#
               #=  if abs(α)<Δ/2
                    h=2/a
                    q=x+2α
                else
                    h=3/a
                    q=(x+h*u)/(1-h*a)
                end =#
                if a*x+u>0
                    h=max(ft-simt,-(a*x+u+a*Δ)/(a*(a*x+u-a*Δ)))# midpoint between asymptote and delta
                  #=   q=x-Δ   #
                    h=h2 =#
                else
                    h=max(ft-simt,-(a*x+u-a*Δ)/(a*(a*x+u+a*Δ)))
                  #=   q=x+Δ 
                    h=h1 =#
                end
                q=(x+h*u)/(1-h*a)#1-h*a non neg because values of h from graph show that h >1/a

            end
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
    #=  if simt>=0.22816661756287676
     dxithr=a*q+u
     @show dxithr
     end =#

    # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
    nextStateTime[i]=simt+h
     return h
end   =#

#= 
function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    
   # a=av[i][i]
    exacteA(qv,cacheA,i,i)
   
    a=cacheA[1]
   quan=quantum[i]
    q=qv[i][0];x=xv[i][0];x1=xv[i][1];
    qaux[i][1]=q
   # olddx[i][1]=x1
    u=x1-a*q
    #uv[i][i][1]=u
    dx=x1
    dxaux[i][1]=x1
    h=0.0
    if a !=0.0
        if dx==0.0
            dx=u+(q)*a
            if dx==0.0
                dx=1e-26
            end
        end
    #for order1 finding h is easy but for higher orders iterations are cheaper than finding exact h using a quadratic,cubic...
    #exacte for order1: h=-2Δ/(u+xa-2aΔ) or h=2Δ/(u+xa+2aΔ)
        h = ft-simt
        q = (x + h * u) /(1 - h * a)
        if (abs(q - x) >  1*quan) # removing this did nothing...check @btime later
          h = (abs( quan / dx));
          q= (x + h * u) /(1 - h * a)
        end
        while (abs(q - x) >  1*quan) 
          h = h * 0.99*(quan / abs(q - x));
          q= (x + h * u) /(1 - h * a)
        end
 
    else
        dx=u
        if dx>0.0
            q=x+quan# 
        else
            q=x-quan
        end
        if dx!=0
        h=(abs(quan/dx))
        else
            h=Inf
        end
    end
    qv[i][0]=q
   # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
   nextStateTime[i]=simt+h
    return h
end   =#

#iterations
#= function updateQ(::Val{2},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{2,Float64}},qaux::Vector{MVector{2,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    exacteA(qv,cacheA,i,i);a=cacheA[1]
   # exacteA(xv,cacheA,i,i);a=cacheA[1]
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
               
        
        if (abs(q - x) > 1.0* quan) # removing this did nothing...check @btime later
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
        while (abs(q - x) >1.0*  quan) && (maxIter>0) && (h>0)
            
          h = h *sqrt(quan / abs(q - x))
        #  h = h *0.99*(1.8*quan / abs(q - x))
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
      #=   if  (abs(q - x) > 2* quan)
            coef=@SVector [quan, -a*quan,(a*a*(x+quan)+a*u1+u2)/2]#
                h1= minPosRoot(coef, Val(2))
              coef=@SVector [-quan, a*quan,(a*a*(x-quan)+a*u1+u2)/2]#
                h2= minPosRoot(coef, Val(2))
  
              if h1<h2
                  h=h1;q=x+quan
              else
                  h=h2;q=x-quan
              end
              qtemp=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
              if abs(qtemp-q)>1e-12
                println("error quad vs qexpression")
              end
              if h1!=Inf && h2!=Inf
                println("quadratic eq double mpr")
              end
  
              if h1==Inf && h2==Inf
                println("quadratic eq NO mpr")
              end


          end  =#
         #= 
                                                                                        if  (abs(q - x) > 2* quan)

                                                                                        # if x2<0.0  #x2 might changed direction...I should use ddx but again ddx=aq+u and q in unknown the sign is unknown
                                                                                              #=   pertQuan=quan-1e-13  #
                                                                                                q=x+pertQuan      # guess q to be thrown up         
                                                                                                coef=@SVector [pertQuan, -a*pertQuan,(a*a*(x+pertQuan)+a*u1+u2)/2]#
                                                                                                h= minPosRoot(coef, Val(2))
                                                                                               pertQuan=quan+1e-13         
                                                                                                coef=@SVector [pertQuan, -a*pertQuan,(a*a*(x+pertQuan)+a*u1+u2)/2]#
                                                                                                h2= minPosRoot(coef, Val(2))
                                                                                                if h2<h
                                                                                                    q=x+pertQuan 
                                                                                                    h=h2
                                                                                                end
                                                                                                ########### q to be thrown down
                                                                                             pertQuan=quan-1e-13
                                                                                                coef=@SVector [-pertQuan, a*pertQuan,(a*a*(x-pertQuan)+a*u1+u2)/2]#
                                                                                                h2= minPosRoot(coef, Val(2))
                                                                                                if h2<h
                                                                                                    q=x-pertQuan 
                                                                                                    h=h2
                                                                                                end
                                                                                                pertQuan=quan+1e-13
                                                                                                coef=@SVector [-pertQuan, a*pertQuan,(a*a*(x-pertQuan)+a*u1+u2)/2]#
                                                                                                h2= minPosRoot(coef, Val(2))
                                                                                                if h2<h
                                                                                                    q=x-pertQuan 
                                                                                                    h=h2
                                                                                                end =#
                                                                                             pertQuan=quan
                                                                                                coef=@SVector [-pertQuan, a*pertQuan,(a*a*(x-pertQuan)+a*u1+u2)/2]#
                                                                                                h= minPosRoot(coef, Val(2))
                                                                                                q=x-pertQuan
                                                                                               #=  if h2<h
                                                                                                    q=x-pertQuan 
                                                                                                    h=h2
                                                                                                end =#
                                                                                                coef=@SVector [pertQuan, -a*pertQuan,(a*a*(x+pertQuan)+a*u1+u2)/2]#
                                                                                                h2= minPosRoot(coef, Val(2))
                                                                                                if h2<h
                                                                                                    q=x+pertQuan 
                                                                                                    h=h2
                                                                                                end
                                                                                        # end 
                                                                                        end =#
                                                                                
             

         

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
                q1=x1
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

    return h
   
end
 =#

#analytic
 function updateQ(::Val{2},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{2,Float64}},qaux::Vector{MVector{2,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    exacteA(qv,cacheA,i,i);a=cacheA[1]
   # exacteA(xv,cacheA,i,i);a=cacheA[1]
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
                                                    
            # quan=quan*1.5
        if a*a*x+a*u1+u2<0.0
            if -(a*a*x+a*u1+u2)/(a*a)<quan # asymptote<delta...no sol ...no need to check
                h=Inf
                q=-(a*u1+u2)/(a*a)
            else
                q=x+quan
                coefi=NTuple{3,Float64}(((a*a*(x+quan)+a*u1+u2)/2,-a*quan,quan))
                h=minPosRootv1(coefi) #needed for dq
            end
        elseif a*a*x+a*u1+u2>0.0
           if -(a*a*x+a*u1+u2)/(a*a)>-quan # asymptote>-delta...no sol ...no need to check
                h=Inf
                q=-(a*u1+u2)/(a*a)
           else
                coefi=NTuple{3,Float64}(((a*a*(x-quan)+a*u1+u2)/2,a*quan,-quan))
                h=minPosRootv1(coefi)
              #  h=mprv2(coefi)
                
                q=x-quan
            end
        else
            q=-u1/a
            h=Inf
        end
         
        if h!=Inf
         α1=1-h*a
        if abs(α1)==0.0
            α1=1e-30*sign(α1)
        end
        q1=(a*q+u1+h*u2)/α1  
        else #h==inf make ddx==0 dq=-u2/a
            q1=-u2/a
        end


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
                q1=x1
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

    return h
   
end  
#= 
function nupdateQ(::Val{2}#= ,cacheA::MVector{1,Float64},map::Function =#,i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},av::Vector{Vector{Float64}},uv::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})where{O}
    #a=getA(Val(Sparsity),cacheA,av,i,i,map)
    a=av[i][i]
    q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent
    olddx[i][1]=x1#appears only here...updated here and used in updateApprox   
    #u1=u1+(simt-tu[i])*u2 # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
    u1=x1-a*qaux[i][1]
    uv[i][i][1]=u1
    uv[i][i][2]=x2-a*q1
    u2=uv[i][i][2]
    #tu[i]=simt  
    # olddx[i][2]=2*x2# 
    ddx=x2  
    quan=quantum[i]
    h=0.0
   #=  if simt == 0.004395600232045285
        @show a
        @show x1
        @show u1
        @show u2

    end =#
    if a!=0.0
        if ddx ==0.0
             ddx=a*a*q+a*u1 +u2
            if ddx==0.0 
                ddx=1e-40# changing -40 to -6 nothing changed
                #println("ddx=0")
            end
        end
        h = ft-simt
        #tempH1=h
        #q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /(1 - h * a + h * h * a * a / 2)
                 q=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
        
        if (abs(q - x) > 2* quan) # removing this did nothing...check @btime later
          h = sqrt(abs(2*quan / ddx)) # sqrt highly recommended...removing it leads to many sim steps..//2* is necessary in 2*quan when using ddx
          q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2)
          
        end
        maxIter=1000
        tempH=h
        while (abs(q - x) > 2* quan) && (maxIter>0) && (h>0)
            
          h = h *sqrt(quan / abs(q - x))
          q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2)
          maxIter-=1
        end

        q1=(a*q+u1+h*u2)/(1-h*a)  #later investigate 1=h*a
    else
        if x2!=0.0  
           h=sqrt(abs(2*quan/x2))   #sqrt necessary with u2
           q=x-h*h*x2/2
           q1=x1+h*x2

        else
           # println("x2==0")
            if x1!=0.0
                h=abs(quan/x1)
                q=x+h*x1
                q1=x1
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

    return h
end =#

function updateQ(::Val{3},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},av::Vector{Vector{Float64}},uv::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})where{O}
    a=av[i][i]
    q=qv[i][0];q1=qv[i][1];q2=2*qv[i][2];x=xv[i][0];x1=xv[i][1];x2=2*xv[i][2];x3=6*xv[i][3];u1=uv[i][i][1];u2=uv[i][i][2];u3=uv[i][i][3]
    elapsed=simt-tq[i]
    qaux[i][1]=q+elapsed*q1+elapsed*elapsed*q2/2#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1+elapsed*q2   ;qaux[i][3]=q2     #never used
    olddx[i][1]=x1  
   # tq[i]=simt
   # elapsed=simt-tu[i]
   # u1=u1+elapsed*u2+elapsed*elapsed*u3/2  
    u1=x1-a*qaux[i][1]
    uv[i][i][1]=u1
   #=  u2=u2+elapsed*u3 
    uv[i][i][2]=u2 =#
   uv[i][i][2]=x2-a*qaux[i][2]  #---------------------------------------------------------------
    u2=uv[i][i][2]
   uv[i][i][3]=x3-a*q2
   u3=uv[i][i][3]
   # tu[i]=simt
    dddx=x3
 
    quan=quantum[i]
    h=0.0
   # println("before q update",abs(q - x) > 2 * quan)
     if a!=0.0
        if dddx ==0.0
            dddx=a*a*a*(q)+a*a*u1+a*u2+u3 #*2
            if dddx==0.0
                dddx=1e-40# changing -40 to -6 nothing changed
                println("dddx=0")  #this appeared once with sys1 liqss3
            end
        end
     

        h = ft-simt
        α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
        λ=x+h*u1+h*h*u2/2+h*h*h*u3/6
        β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+λ
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
        q = β/γ
        if (abs(q - x) >  2*quan) # removing this did nothing...check @btime later
          h = cbrt(abs((6*quan) / dddx));
          #h= cbrt(abs((q-x) / x3));#h=cbrt(abs(6*(q-x) / x3))# shifts up a little
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
        q = β/γ
        end

        
        maxIter=515
        while (abs(q - x) >  2*quan) && (maxIter>0)
            maxIter-=1
         # h = h *(0.98*quan / abs(q - x));
          h = h *cbrt(quan / abs(q - x));
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)

        q = β/γ
        end
        q1=(a*(1-h*a)*q+u1*(1-h*a)-h*h*(a*u2+u3)/2)/(1-h*a+h*h*a*a/2)
        q2=(a*q1+u2+h*u3)/(1-h*a)

        if maxIter <200
            @show maxIter
        end
        if h==0.0
            @show h
        end
 
    else
       
        if x3!=0.0
            h=cbrt(abs(6*quan/x3))
            q=x+h*h*h*x3/6
            q1=x1-x3*h*h/2   #*2
            q2=x2+h*x3
       else
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
                    h=Inf
                    
                end
                q1=x1
                
            end
            q2=x2
        end 
    end
    qv[i][0]=q
    qv[i][1]=q1 
    qv[i][2]=q2/2  
    nextStateTime[i]=simt+h
    return nothing
end

function nupdateQ(::Val{3},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},av::Vector{Vector{Float64}},uv::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})where{O}
    a=av[i][i]
    q=qv[i][0];q1=qv[i][1];q2=2*qv[i][2];x=xv[i][0];x1=xv[i][1];x2=2*xv[i][2];x3=6*xv[i][3];u1=uv[i][i][1];u2=uv[i][i][2];u3=uv[i][i][3]
    elapsed=simt-tq[i]
    qaux[i][1]=q+elapsed*q1+elapsed*elapsed*q2/2#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1+elapsed*q2   ;qaux[i][3]=q2     #never used
    olddx[i][1]=x1  
   # tq[i]=simt
    #elapsed=simt-tu[i]
    #u1=u1+elapsed*u2+elapsed*elapsed*u3/2  
    u1=x1-a*qaux[i][1]
    uv[i][i][1]=u1
   # u2=u2+elapsed*u3 
   # uv[i][i][2]=u2
    uv[i][i][2]=x2-a*qaux[i][2]
    u2=uv[i][i][2]
   uv[i][i][3]=x3-a*q2
   u3=uv[i][i][3]
    #tu[i]=simt
    dddx=x3
   
    quan=quantum[i]
    h=0.0

     if a!=0.0
        if dddx ==0.0
            dddx=a*a*a*(q)+a*a*u1+a*u2+u3 #*2
            if dddx==0.0
                dddx=1e-40# changing -40 to -6 nothing changed
               println("nupdate dddx=0")  #this appeared once with sys1 liqss3
            end
        end
 
        h = ft-simt
        α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
        λ=x+h*u1+h*h*u2/2+h*h*h*u3/6
        β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+λ
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
        q = β/γ
        if (abs(q - x) > 2* quan) # removing this did nothing...check @btime later
          h = cbrt(abs((6*quan) / dddx));
         
          #h= cbrt(abs((q-x) / x3));#h=cbrt(abs(6*(q-x) / x3))# shifts up a little
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
        q = β/γ
        end

        
        maxIter=515
        while (abs(q - x) > 2* quan) && (maxIter>0)
            maxIter-=1
         # h = h *(0.98*quan / abs(q - x));
          h = h *cbrt(quan / abs(q - x))
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)

        q = β/γ
        end
        q1=(a*(1-h*a)*q+u1*(1-h*a)-h*h*(a*u2+u3)/2)/(1-h*a+h*h*a*a/2)
        q2=(a*q1+u2+h*u3)/(1-h*a)
        if maxIter <200
            @show maxIter
        end
        if h==0.0
            @show h
        end
 
    else
       
        if x3!=0.0
            h=cbrt(abs(6*quan/x3))
            q=x+h*h*h*x3/6
            q1=x1-x3*h*h/2   #*2
            q2=x2+h*x3
       else
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
        end 
    end
    qv[i][0]=q
    qv[i][1]=q1 
    qv[i][2]=q2/2  
    nextStateTime[i]=simt+h
    return nothing
end


#= function Liqss_ComputeNextTime(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64})where{T}
    #= q=qv[i][0];x=xv[i][0];q1=qv[i][1];x1=xv[i][1];x2=xv[i][2]
    coef=@SVector [q - x , q1-x1,-x2] =#
    coef=@SVector [qv[i][0]- xv[i][0] , -xv[i][1],]#
    nextStateTime[i] = simt + minPosRoot(coef, Val(1))
end =#
#= function Liqss_ComputeNextTime1(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64})where{T}
    q=qv[i][0];x=xv[i][0];x1=xv[i][1]
    #if xv[i][1] !=0.0
    if  x1!=0.0  
        nextStateTime[i]=simt+(abs((q-x)/(x1)))  
    else
        nextStateTime[i]=Inf
    end
end =#
function Liqss_ComputeNextTime(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64})
    q=qv[i][0];x=xv[i][0];x1=xv[i][1]
    if  x1!=0.0  
        nextStateTime[i]=simt+(abs((quantum[i])/(x1)))  
    else
        nextStateTime[i]=Inf
    end
end

function Liqss_ComputeNextTime(::Val{2}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64})
    q=qv[i][0];x=xv[i][0];q1=qv[i][1];x1=xv[i][1];x2=xv[i][2]
    if  x2!=0.0  
        nextStateTime[i]=simt+sqrt(abs((q-x)/(x2)))  
    else
        nextStateTime[i]=Inf
    end
end

function Liqss_reComputeNextTime(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,a::Vector{Vector{Float64}} =#)
    dt=0.0; q=qv[i][0];x=xv[i][0];x1=xv[i][1]
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
    if nextStateTime[i]<simt # this is coming from the fact that a variable can reach 2quan distance when it is not its turn, then computation above gives next=simt+(p-p)/dx...p-p should be zero but it can be very small negative
        nextStateTime[i]=simt+1e-12
    end
end



function Liqss_reComputeNextTime(::Val{2}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,a::Vector{Vector{Float64}} =#)
    q=qv[i][0];x=xv[i][0];q1=qv[i][1];x1=xv[i][1];x2=xv[i][2];quani=quantum[i]
    β=0
    coef=@SVector [q-x#= -1e-13 =#, q1-x1,-x2]#

           #=  nextStateTime[i]=simt + minPosRoot(coef, Val(2))
            #= coef=setindex(coef, q-x-1e-13,1)
            timetemp = simt + minPosRoot(coef, Val(2))
            if timetemp < nextStateTime[i] 
                nextStateTime[i]=timetemp
            end =#
            coef=setindex(coef, q-x+1e-13,1)
            timetemp = simt + minPosRoot(coef, Val(2))
            if timetemp < nextStateTime[i] 
                nextStateTime[i]=timetemp
            end
            coef=setindex(coef, q-x,1) =#
            nextStateTime[i] = simt + minPosRoot(coef, Val(2))
           #=  if timetemp < nextStateTime[i] 
                nextStateTime[i]=timetemp
            end =#
          #=   h=sqrt(abs((q-x)/x2))#2delta/2 x2
            α=q-x+(q1-x1)*h-x2*h*h
            maxIter=500
            while abs(α)>1e-13 && maxIter>0
                maxIter-=1
                h=h*0.96*1e-13/abs(α)
                α=q-x+(q1-x1)*h-x2*h*h
            end
            if maxIter>0
                nextStateTime[i]=simt+h
            else
                nextStateTime[i]=Inf
            end =#
        if q-x >0.0#1e-9
            coef=setindex(coef, q-x-2*quantum[i],1)
            timetemp = simt + minPosRoot(coef, Val(2))
            if timetemp < nextStateTime[i] 
                nextStateTime[i]=timetemp
            end
        elseif  q-x <0.0#-1e-9
            coef=setindex(coef, q-x+2*quantum[i],1)
            timetemp = simt + minPosRoot(coef, Val(2))
            if timetemp < nextStateTime[i] 
                nextStateTime[i]=timetemp
            end
        else
          #  nextStateTime[i]=simt+Inf#1e-19
            #nextStateTime[i]=simt+1e-19
          #=  if q-x==0.0
            nextStateTime[i]=simt+Inf#1e-19 #
           else
            nextStateTime[i]=simt+1e-12#Inf#1e-19 #
           end =#
        end


        if nextStateTime[i]<simt # this is coming from the fact that a variable can reach 2quan distance when it is not its turn, then computation above gives next=simt+(p-p)/dx...p-p should be zero but it can be very small negative
            nextStateTime[i]=simt+1e-15
           # @show simt,nextStateTime[i],i,x,q,quantum[i],xv[i][1]
        end


end
 



function Liqss_reComputeNextTime(::Val{3}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,a::Vector{Vector{Float64}} =#)
    q=qv[i][0];x=xv[i][0];q1=qv[i][1];x1=xv[i][1];x2=xv[i][2];q2=qv[i][2];x3=xv[i][3]
    coef=@SVector [q - x , q1-x1,q2-x2,-x3]# x and q might get away even further(change of sign) and we want that to be no more than another quan
    nextStateTime[i] = simt + minPosRoot(coef, Val(3))
   
    if q-x >0
        coef=setindex(coef, q-x-2*quantum[i],1)
        timetemp = simt + minPosRoot(coef, Val(3))
        if timetemp < nextStateTime[i] 
            nextStateTime[i]=timetemp
        end
    elseif  q-x <0
        coef=setindex(coef, q-x+2*quantum[i],1)
        timetemp = simt + minPosRoot(coef, Val(3))
        if timetemp < nextStateTime[i] 
            nextStateTime[i]=timetemp
        end
    else
        nextStateTime[i]=simt+Inf#1e-19
    end

    
end


#######################################################################################################################################################
#= function updateLinearApprox(::Val{1},i::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},simt::Float64)where{T,O}
    diffQ=q[i][0]-qaux[i][1]
    if diffQ != 0.0
        a[i][i]=(x[i][1]-olddx[i][1])/diffQ
    else
        a[i][i]=0.0
    end
    u[i][i][1]=x[i][1]-a[i][i]*q[i][0]   #if a==0 u same as derx meaning that in updateQ if derx> we dont have to check if u>0 ....note1
    return nothing
end =#

function updateLinearApprox(i::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}})where{O}
    diffQ=q[i][0]-qaux[i][1]
     if diffQ != 0.0
        a[i][i]=(x[i][1]-olddx[i][1])/diffQ
    else
        a[i][i]=0.0
    end
    return nothing
end
#nupdateLinear differ from updateLinear in that nupdate uses with qminus instead of qaux=qminus+e*dq
function nupdateLinearApprox(i::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},qminus::Vector{Float64},olddx::Vector{MVector{O,Float64}})where{O}
    diffQ=q[i][0]-qminus[i]
      if diffQ != 0.0
       a[i][i]=(x[i][1]-olddx[i][1])/diffQ
    else
        a[i][i]=0.0
    end
    return nothing
end


function updateOtherApprox(k::Int,j::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}})where{O}
    diffQ=q[j][0]-qaux[j][1]
    if diffQ != 0.0
      a[k][j]=(x[k][1]-olddx[k][1])/diffQ
    else
      a[k][j]=0.0
    end
    return nothing
  end



  #= function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    
    # a=av[i][i]
     exacteA(qv,cacheA,i,i)
    
     a=cacheA[1]
 

  #=    if i==1
        if abs(a+1.1)>1e-3
            @show i,a
        end
    else
        if abs(a+20.0)>1e-3
            @show i,a
        end
    end =#

     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
    #=  olddx[i][1]=x1 =#
     u=x1-a*q
     #uv[i][i][1]=u
     dx=x1
     dxaux[i][1]=x1
     h=0.0
     if a !=0.0
         if dx==0.0
             dx=u+(q)*a
             if dx==0.0
                 dx=1e-26
             end
         end
     #for order1 finding h is easy but for higher orders iterations are cheaper than finding exact h using a quadratic,cubic...
     #exacte for order1: h=-2Δ/(u+xa-2aΔ) or h=2Δ/(u+xa+2aΔ)
        #=  h = (ft+100.0-simt)
         q = (x + h * u) /(1 - h * a)
         if (abs(q - x) >  quantum[i]) =#
            h1=-1.0
            h0 = (ft-simt)
            q0 = (x + h0 * u) /(1 - h0 * a)
            if (abs(q - x) <  quantum[i])
                h1=h0
            end
        # end
        #=  if (abs(q - x) >  2*quantum[i]) # removing this did nothing...check @btime later
           h = (abs( quantum[i] / dx));
           q= (x + h * u) /(1 - h * a)
         end =#
        #=  while (abs(q - x) >  2*quantum[i]) 
           h = h * 1.98*(quantum[i] / abs(q - x));
           q= (x + h * u) /(1 - h * a)
         end =#
         coefQ=1
         #if (abs(q - x) >  coefQ*quantum[i])
            h2=coefQ*quantum[i]/(a*(x+coefQ*quantum[i])+u)
            
            if h2<0
                h2=-coefQ*quantum[i]/(a*(x-coefQ*quantum[i])+u)
                q=x-coefQ*quantum[i]
            else
                q=x+coefQ*quantum[i]
            end
           
         #end
        if h2<0
            h=Inf
        end
        h=max(h1,h2)
       #=   if h==ft+100.0-simt
            @show i,simt,x,x1,q,a*q+u
            
         end =#

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
=#
 
    