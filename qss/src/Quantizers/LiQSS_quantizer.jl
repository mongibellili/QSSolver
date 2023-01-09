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
function updateQ(::Val{1},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64, nextTime::MVector{T,Float64})where{T,O}
    qaux[i][1]=qv[i][0]# index shift....sorry but be careful: taylor 1st elemtn is at 0, a vect 1st elemnt is at 1
    olddx[i][1]=xv[i][1]
    #q[i][0]=x[i][0]
    q=qv[i][0]
    u=uv[i][i][1]  # for order 2: u=u+tu*deru
    #dq=0.0
    x=xv[i][0]
    dx=xv[i][1]
    a=av[i][i]
    h=0.0
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
        if (abs(q - x) >  quantum[i]) # removing this did nothing...check @btime later
          h = (abs( quantum[i] / dx));
          q= (x + h * u) /(1 - h * a)
        end
        while (abs(q - x) >  quantum[i]) 
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
   nextTime[i]=simt+h
    return nothing
end
function updateQ(::Val{2},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64, nextTime::MVector{T,Float64})where{T,O}
    q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
   

    #= if simt>400
        @show simt
       end =#

  # println("qaux i inside updateQ= ",qaux[i][1])
    #q=qaux[i][1]# not needed...q used only in approx ddx which not needed to be exacte
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent
   # tq[i]=simt
    olddx[i][1]=x1#appears only here...updated here and used in updateApprox   
    u1=u1+(simt-tu[i])*u2 # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
    uv[i][i][1]=u1
    tu[i]=simt  
     olddx[i][2]=2*x2# 
    ddx=x2
    a=av[i][i]

   #=  if abs(av[1][1])>1e-3 || abs(av[1][2]-0.01)>1e-3 || abs(av[2][1]+100)>1e-3 || abs(av[2][2]+100)>1e-3  #sys2 testing
        @show simt
        @show i
        @show av
      end =#

      #= if abs(uv[1][2][1])>1e-3 || abs(uv[2][1][1]-2020)>1e-3 
        @show simt
        @show i
        @show uv[1][1][1],xv[1][1]
        @show uv[1][2][1]
        @show uv[2][1][1]
        @show uv[2][2][1],xv[2][1]
      end =#

    quan=quantum[i]
    h=0.0
    if a!=0.0
        if ddx ==0.0
            #= ddx=a*a*q+a*u1 +u2
            if ddx==0.0 =#
                ddx=1e-26# changing -40 to -6 nothing changed
                println("ddx=0")
            #end
        end
        h = ft-simt
        #tempH1=h
        q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                 (1 - h * a + h * h * a * a / 2)
<<<<<<< HEAD
        
        if (abs(q - x) >  2*quan) # removing this did nothing...check @btime later
=======
        if (abs(q - x) >  quan) # removing this did nothing...check @btime later
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
          h = sqrt(abs(2*quan / ddx)) # sqrt highly recommended...removing it leads to many sim steps..//2* is necessary in 2*quan when using ddx
          q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2)
           #= if h>ft-simt
            @show h
           end =#
        end
<<<<<<< HEAD
        maxIter=1000
        tempH=h
        while (abs(q - x) >  2*quan) && (maxIter>0)
            
          h = h *1.96*(quan / abs(q - x))
=======
        while (abs(q - x) >  quan) 
          h = h *quan / abs(q - x)
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
          q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2)
          maxIter-=1
         #=  if maxIter < 990
                    println("maxiter of updateQ      = ",maxIter)
                    @show  h, q-x
                   
           end =#
        end
        if maxIter < 900
             println("maxiter of updateQ val2      = ",maxIter)
            @show tempH, h, simt
        end
       #=  @show tempH1,tempH2
        @show q-x
        @show i,h, a, u1, u2  =#
        q1=(a*q+u1+h*u2)/(1-h*a)  #later investigate 1=h*a
       #=  @show q1
        
        @show -1.5775988642611344e-20-1.0842021724855044e-19+q-x+h*(q1-(a*q+u1))-h*h*(a*q1+u2)/2
        coef=@SVector [-1.5775988642611344e-20-1.0842021724855044e-19+q-x , q1-(a*q+u1),-(a*q1+u2)/2]#
      
        hi = 0.0 + minPosRoot(coef, Val(2))
        @show hi =#
        if abs(q-x)>2*quantum[i]#uncommenting this did nothing
            println("at end of q update if side q-x >quan !!!! simt= ",simt)
        end

    else
       # println("a=0")
        #ddx=u2
        #= if x2>0.0  #if ddx>0.0 same results
            q=x+quan
        else# elseif x2<0   ......else q=x??
            q=x-quan
        end =#
        #q=x+quantum[i]  #works fine !!!
       # q=x-quantum[i] #errors ...solution escapes up...later test if this behavior is specific to this problem or it is a general thing
       # q=x+2*quantum[i]  #2*Δ errors solution escapes down
       # q=x+quantum[i]   #removing it errors
        if x2!=0.0
          
            h=sqrt(abs(2*quan/x2))   #sqrt necessary with u2
           # q=x+h*x1+h*h*x2/2
           q=x-h*h*x2/2
           q1=x1+h*x2
          # q=x-h*h*u2/2
           # q1=u1+h*u2  #(250 allocations: 16.07 KiB)
          #  @show h,x2,i,u1,u2
        else
            println("x2==0")
            if x1!=0.0
                h=abs(quan/x1)
                q=x+h*x1
                q1=x1
            else
                h=Inf
                q=x
                q1=x1
            end
        #=q1=x1+h*x2  #(252 allocations: 16.07 KiB)
        else
            q1=x1 =#
            #println("ddx=0")
        end 
        if abs(q-x)>2*quantum[i]+1e-12
            println("at end of q update: else side, q-x >quan !!!! simt= ",simt,", q-x= ",q-x,", quan= ",quantum[i])
        end
    end
    
    #olddx[i][2]=ddx  #olddx[i][2] never used again so no need to update it 
    qv[i][0]=q
    qv[i][1]=q1  
   # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
   nextTime[i]=simt+h
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
#= function updateQ(::Val{2},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    q=xv[i][0]
    q1=qv[i][1]
    x=xv[i][0]
    x1=xv[i][1]
    x2=xv[i][2]
    u=uv[i][i][1]
    u1=uv[i][i][2]
    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent
    tq[i]=simt
    olddx[i][1]=x1#appears only here...updated here and used in updateApprox   
    u=u+(simt-tu[i])*u1 # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
    uv[i][i][1]=u
    tu[i]=simt
    
   # olddx[i][2]=2*x2# 
    ddx=x2*2
    dq=0.0
    a=av[i][i]
    if a!=0.0
        if x2>0
        #if ddx ==0.0
            ddx=a*a*(q+quantum[i])+a*u +u1
            if ddx>=0.0
                dq=quantum[i]# changing -40 to -6 nothing changed
            else
                dq=-u1/(a*a)-u/(a)-q
            end
            if abs(dq)>quantum[i]
                dq=quantum[i]
            end
        else
            ddx=a*a*(q-quantum[i])+a*u +u1
            if ddx<=0.0
                dq=-quantum[i]# changing -40 to -6 nothing changed
            else
                dq=-u1/(a*a)-u/(a)-q
            end
            if abs(dq)>quantum[i]
                dq=-quantum[i]
            end
        end
        
        q1=a*q+u  # this is missing +h*x2=+h(a*q1+u1)...as if there is an order drop here...for this version stating q1=a*oldq+u which is old x1 ...trying to make it future leads to nextline error
        q+=dq # in this approach new q has to be after q1...ie q1 cannot be calculated using new q....infinite loop..program does not terminate
 
    else
        
        if x2>0.0
            dq=-quantum[i]
        else
            dq=quantum[i]
        end
        q+=dq
       # q=x+quantum[i]  #works fine !!!
       # q=x-quantum[i] #errors ...solution escapes up...later test if this behavior is specific to this problem or it is a general thing
       # q=x+2*quantum[i]  #2*Δ errors solution escapes down
       # q=x+quantum[i]   #removing it errors
       if x2!=0.0
            h=sqrt(abs(2*quantum[i]/(x2)))
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
    olddx[i][2]=ddx  #olddx[i][2] never used again so no need to update it 
    qv[i][0]=q
    qv[i][1]=q1  
    return nothing
end  =#

#= function updateQ(::Val{2},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    q=xv[i][0]  # q=x used in liqss1 and liqss3
    q1=qv[i][1]
   # q2=qv[i][2]
    #----------------
    x=xv[i][0]
    x1=xv[i][1]
    x2=xv[i][2]*2
   # x3=xv[i][3]*6
    #----------------
    u=uv[i][i][1]
    u1=uv[i][i][2]
   # u2=uv[i][i][3]
    #----------------
    qaux[i][1]=q+(simt-tq[i])*q1#+(simt-tq[i])*(simt-tq[i])*q2#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1#+(simt-tq[i])*q2                     #never used
    #qaux[i][3]=q2                                     #never used
    tq[i]=simt
    olddx[i][1]=x1 #olddx[i][1] appears only here...updated here and used in updateApprox
    
    u=u+(simt-tu[i])*u1#+(simt-tu[i])*(simt-tu[i])*u2/2  # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
    uv[i][i][1]=u
   # u1=u+(simt-tu[i])*u2 #2*u2
   # uv[i][i][2]=u1
    tu[i]=simt
    
   # olddx[i][2]=2*x2#  olddx[i][2] not used elsewhere....so not needed
   # ddx=olddx[i][2]
   dq=0.0
    a=av[i][i]
    if a!=0.0
        if x2>0
            dddx=a*a*(q+quantum[i])+a*u+u1 #*2
            if dddx >=0 
                dq=quantum[i]
            else
                dq=-u1/(a*a)-u/a-q  #*2
            end
            if abs(dq)>quantum[i]
                dq=quantum[i]
            end
        else
            dddx=a*a*(q-quantum[i])+a*u+u1  #*2
            if dddx <=0 
                dq=-quantum[i]
            else
                dq=-u1/(a*a)-u/a-q  #*2
            end
            if abs(dq)>quantum[i]
                dq=-quantum[i]
            end
        end
        q1=a*q+u
       # q2=a*qv[i][1]+u1
        #= if q1*x1 <0.0 && !flag3
            dq=qaux[i][1]-q-abs(quantum[i])*0.1
        else
            dq=qaux[i][1]-q+abs(quantum[i])*0.1
        end =#
    else #a==0
        if x2>0.0
            dq=-quantum[i]
        else
            dq=+quantum[i]
        end
        if x2!=0.0
            h=sqrt(abs(4*quantum[i]/(x2)))
            q1=u+h*u1#+h*h*u2/2   #*2
           # q2=u1+h*u2  #*2
       else
           #=  q1=x1#u
            q2=x2#u1 =#
            q1=u
            #q2=u1
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
   # qv[i][2]=q2/2   #/2???

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

function updateQ(::Val{3},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64, nextTime::MVector{T,Float64})where{T,O}
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
    h=0.0
   # println("before q update",abs(q - x) > 2 * quan)
     if a!=0.0
        if dddx ==0.0
           # dddx=a*a*a*(q)+a*a*u+a*u1+u2 #*2
           # if dddx==0.0
                dddx=1e-26# changing -40 to -6 nothing changed
<<<<<<< HEAD
                println("dddx=0")  #this appeared once with sys1 liqss3
=======
                println("dddx=0")
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
           # end
        end
        h = ft-simt
        α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
        β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
<<<<<<< HEAD
      #=  β=h*h*(2*a*u1+h*a*(2*u2/3-4*a*u1/3+u1/2+h*(a*a*u1/3+u3/6-11*a*u2/12+h*a*(a*u2/4-u3/4+h*a*u3/12))))/(1-a*h)+x+h*h*h*u3/6
       γ=1-a*h+h*h*a*a/2-h*h*h*a*a*a/6 =#
        q = β/γ
        if (abs(q - x) >  quan) # removing this did nothing...check @btime later
          h = cbrt(abs((6*quan) / dddx));
          #h= cbrt(abs((q-x) / x3));#h=cbrt(abs(6*(q-x) / x3))# shifts up a little
          
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        #γ=(1-a*h+h*h*a*a/2-h*h*h*a*a*a/6)/(1-a*h+h*h*a*a/2)
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)


        #= β=h*h*(2*a*u1+h*a*(2*u2/3-4*a*u1/3+u1/2+h*(a*a*u1/3+u3/6-11*a*u2/12+h*a*(a*u2/4-u3/4+h*a*u3/12))))/(1-a*h)+x+h*h*h*u3/6
       γ=1-a*h+h*h*a*a/2-h*h*h*a*a*a/6 =#
=======
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
        q = β/γ
        if (abs(q - x) >  quan) # removing this did nothing...check @btime later
            h = cbrt(abs((6*quan) / dddx));
            α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
            β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
            γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
            q = β/γ
        end
<<<<<<< HEAD
       # println(abs(q - x) > 2 * quan)
        
        maxIter=515
        while (abs(q - x) >  quan) && (maxIter>0)
            maxIter-=1
          h = h *(0.95*quan / abs(q - x));
=======
        maxIter=215
        while (abs(q - x) >  quan) && (maxIter>0)
          maxIter-=1
          h = h *quan / abs(q - x);
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
          γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
          q = β/γ
        end
         if maxIter < 200
             println("maxiter of mpudate= ",maxIter)
<<<<<<< HEAD
         end


=======
         end =#
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
        q1=(a*(1-h*a)*q+u1*(1-h*a)-h*h*(a*u2+u3)/2)/(1-h*a+h*h*a*a/2)
        q2=(a*q1+u2+h*u3)/(1-h*a)
    else
        if x3!=0.0
            h=cbrt(abs(6*quan/x3))
            q=x+h*h*h*x3/6
            q1=x1-x3*h*h/2   #*2
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
    nextTime[i]=simt+h
    return nothing
end

<<<<<<< HEAD
#= function updateQ(::Val{3},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    q=qv[i][0];q1=qv[i][1];q2=2*qv[i][2];x=xv[i][0];x1=xv[i][1];x2=2*xv[i][2];x3=6*xv[i][3];u1=uv[i][i][1];u2=uv[i][i][2];u3=uv[i][i][3]
    qaux[i][1]=q+(simt-tq[i])*q1+(simt-tq[i])*(simt-tq[i])*q2/2#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1+(simt-tq[i])*q2   ;qaux[i][3]=q2     #never used
    olddx[i][1]=x1#+(simt-tq[i])*x2 #later change to tx....appears only here...updated here and used in updateApprox   
    tq[i]=simt
    u1=u1+(simt-tu[i])*u2+(simt-tu[i])*(simt-tu[i])*u3/2  
    uv[i][i][1]=u1
    u2=u2+(simt-tu[i])*u3 #2*u2
    uv[i][i][2]=u2
    tu[i]=simt
    
   # olddx[i][2]=2*x2# 
    dddx=x3
    a=av[i][i]
    quan=quantum[i]
     if a!=0.0
        if dddx ==0.0
           # dddx=a*a*a*(q)+a*a*u+a*u1+u2 #*2
           # if dddx==0.0
                dddx=1e-26# changing -40 to -6 nothing changed
           # end
        end
         #=   coef=@SVector [ -2* quantum[i], x*a+2*quantum[i]+u,(a*u+u1+x*a*a+2*quantum[i]*a*a)/2]#*2
        h =  minPosRoot(coef, Val(2))
        if h==Inf
             coef=@SVector [ 2* quantum[i], x*a-2*quantum[i]+u,(a*u+u1+x*a*a-2*quantum[i]*a*a)/2]#*2
             h =  minPosRoot(coef, Val(2))
        end
        q=(x+h*u+h*h*(a*u+u1)/2)/(1-h*a-h*h*a*a/2)
        q1=(a*q+u+h*u1)/(1-h*a) =#
        #=         coef=@SVector [ -2* quantum[i], x*a+4*a*quantum[i]+u,(-a*u+u1-x*a*a-2*quantum[i]*a*a)/2]#*2
        h =  minPosRoot(coef, Val(2))
        if h==Inf
             coef=@SVector [ 2* quantum[i], x*a-4*a*quantum[i]+u,(-a*u+u1-x*a*a+2*quantum[i]*a*a)/2]#*2
             h =  minPosRoot(coef, Val(2))
        end
        q=((x+h*u+h*h*u1/2)*(1-h*a)+(u+h*u1)*h*h*a/2)/(1-2*h*a+h*h*a*a/2) =#

        h = ft-simt
        q = ((x + h * u1 + h * h / 2 * (a*u1+u2))  + (h*h * h / 6 ) * (a*a*u1 + a*u2+u3)) /
                 (1 - h * a + h * h * a * a / 2-h*h*h*a*a*a/6)
        if (abs(q - x) > 2 * quan) # removing this did nothing...check @btime later
          h = cbrt(abs((quan) / dddx));
          #h= cbrt(abs((q-x) / x3));#h=cbrt(abs(6*(q-x) / x3))# shifts up a little
          
          q = ((x + h * u1 + h * h / 2 * (a*u1+u2))  + (h*h * h / 6 ) * (a*a*u1 + a*u2+u3)) /
                 (1 - h * a + h * h * a * a / 2-h*h*h*a*a*a/6)
        end
        
        while (abs(q - x) > 2 * quan) 
          h = h *sqrt(quan / abs(q - x));
          q = ((x + h * u1 + h * h / 2 * (a*u1+u2))  + (h*h * h / 6 ) * (a*a*u1 + a*u2+u3)) /
                 (1 - h * a + h * h * a * a / 2 -h*h*h*a*a*a/6)
        end


        q1=(a*q+u1+h*u2+h*h*(a*u2+u3)/2)/(1-h*a-h*h*a*a/2)
        q2=(a*q1+u2+h*u3)/(1-h*a)
 
    else
        #dddx=u2
        if x3>0.0
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
       end 
    end
    #= if abs(q-x)>2*quantum[i]#uncommenting this did nothing
        q=x
    end =#
    #olddx[i][2]=ddx  #olddx[i][2] never used again so no need to update it 
    qv[i][0]=q
    qv[i][1]=q1 
    qv[i][2]=q2/2  
    return nothing
end =#



#= function Liqss_ComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    #= q=qv[i][0];x=xv[i][0];q1=qv[i][1];x1=xv[i][1];x2=xv[i][2]
    coef=@SVector [q - x , q1-x1,-x2] =#
   #= if i==1
    @show qv[i][0]- xv[i][0]#-1.000020000008206e-5
    @show qv[i][1]- xv[i][1]# 0.01819692125101185
   end =#
    coef=@SVector [qv[i][0]- xv[i][0] , -xv[i][1],]#
    nextTime[i] = currentTime + minPosRoot(coef, Val(1))
end =#
#= function Liqss_ComputeNextTime1(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    q=qv[i][0]
    x=xv[i][0]
   
    x1=xv[i][1]
   
    #if xv[i][1] !=0.0
    if  x1!=0.0  
        nextTime[i]=currentTime+(abs((q-x)/(x1)))  
    else
        nextTime[i]=Inf
    end
end =#
function Liqss_ComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    q=qv[i][0]
    x=xv[i][0]
   
    x1=xv[i][1]
   
    #if xv[i][1] !=0.0
    if  x1!=0.0  
        nextTime[i]=currentTime+(abs((quantum[i])/(x1)))  
    else
        nextTime[i]=Inf
    end
end

#= function Liqss_ComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    #= q=qv[i][0];x=xv[i][0];q1=qv[i][1];x1=xv[i][1];x2=xv[i][2]
    coef=@SVector [q - x , q1-x1,-x2] =#
   #= if i==1
    @show qv[i][0]- xv[i][0]#-1.000020000008206e-5
    @show qv[i][1]- xv[i][1]# 0.01819692125101185
   end =#
    coef=@SVector [qv[i][0]- xv[i][0] , qv[i][1]-xv[i][1],-xv[i][2]]#
    nextTime[i] = currentTime + minPosRoot(coef, Val(2))
end =#
=======

>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
function Liqss_ComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
<<<<<<< HEAD
    #if xv[i][1] !=0.0
    if  x2!=0.0  
        nextTime[i]=currentTime+sqrt(abs((q-x)/(x2)))  
    else
        nextTime[i]=Inf
    end
end
#= function Liqss_ComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
=======
    
   
    coef=@SVector [q - x , q1-x1,-x2]#
    nextTime[i] = currentTime + minPosRoot(coef, Val(2))
    
    #nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end
function Liqss_ComputeNextTime(::Val{3}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64})where{T}
>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
<<<<<<< HEAD
    #if xv[i][1] !=0.0
    if  x2!=0.0  
        nextTime[i]=currentTime+sqrt(abs((quantum[i])/(x2)))  
    else
        nextTime[i]=Inf
    end
end  =#
=======
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




>>>>>>> ee1ea40e28192d3d5a4e88a58aa57a201afb63b3

##########################################################################################################################################################
#= function Liqss_reComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
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
end =#
function Liqss_reComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    dt=0.0
    q=qv[i][0]
    x=xv[i][0]
    if xv[i][1] !=0.0 #&& abs(q-x)>quantum[i]/10
        dt=(q-x)/xv[i][1]
        if dt>0.0
            nextTime[i]=currentTime+dt# later guard against very small dt
        elseif dt<0.0
            if xv[i][1]>0.0  
                nextTime[i]=currentTime+(quantum[i])/xv[i][1]
            else
                nextTime[i]=currentTime+(-quantum[i])/xv[i][1]
            end
           #=  if nextTime[i] < currentTime 
                println("liqss1 recompute: x heading neither towards nor 2delta away! at time= ",currentTime) 
                nextTime[i]=currentTime+1e-1
            end =#
        else#q=x equilibruim point reached
            nextTime[i]=Inf# for simple recompute maube nexttime=now
            println("dt==0")
        end
    else
        nextTime[i]=Inf
    end

end


#= function Liqss_reComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
    #if xv[i][1] !=0.0
    if  x2!=0.0 #&& a[i][i] != 0.0  # a!=0 cuz when a==0 x1=0+u  and then u updated below as u=x1-0 meaning that x1==u=oldx1....no change ie nexttime should be Inf???
       # nextTime[i]=currentTime+sqrt(abs((q-x)/(x2)))  #64.829 μs (580 allocations: 46.54 KiB)
     #  nextTime[i]=currentTime+sqrt(abs((quantum[i])/(x2)))  #64.503 μs (584 allocations: 46.66 KiB)
       #nextTime[i]=currentTime+(abs((x1)/(x2)))#nextTime[i]=currentTime+sqrt(abs((x1)/(x2))) # 54.390 μs (584 allocations: 46.66 KiB)
      # nextTime[i]=currentTime+(abs((q1-x1)/(x2)))# 57.187 μs (586 allocations: 46.73 KiB)
       coef=@SVector [#= q - x =# quantum[i] , q1-x1,-x2]#135.287 μs (580 allocations: 46.54 KiB)
       nextTime[i] = currentTime + minPosRoot(coef, Val(2)) 
    else
        nextTime[i]=Inf
    end
    #removing the previous 6 lines and adding time3 -->78.725 μs (656 allocations: 48.91 KiB)
    coef=@SVector [q-x+2*quantum[i], q1-x1,-x2]#
    time1 = currentTime + minPosRoot(coef, Val(2))
    coef=setindex(coef, q-x-2* quantum[i],1)
    time2 = currentTime + minPosRoot(coef, Val(2))
    time1 = time1 < time2 ? time1 : time2    
   #=coef=setindex(coef,q - x,1)
    time3 = currentTime + minPosRoot(coef, Val(2))
    time1 = time1 < time3 ? time1 : time3  =#  
    nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    #= if q*q1<0 && abs(q) > 10.0*quantum[i] # uncomment did nothing
        time3=-(q/q1)-2*abs(quantum[i]/q1)
        if time3>0
        nextTime[i] = currentTime+time3 < nextTime[i] ? currentTime+time3 : nextTime[i]
        end
    end  =#   
   #=  else
        nextTime[i]=Inf
    end =#





end =#
function Liqss_reComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
    if 350.95<currentTime<350.96
        @show i,x,q,currentTime
    end
    coef=@SVector [q-x, q1-x1,-x2]#
        nextTime[i]=currentTime + minPosRoot(coef, Val(2))
        if q-x >0
            coef=setindex(coef, q-x-2*quantum[i],1)
            timetemp = currentTime + minPosRoot(coef, Val(2))
            if timetemp < nextTime[i] 
                nextTime[i]=timetemp
            end
        elseif  q-x <0
            coef=setindex(coef, q-x+2*quantum[i],1)
            timetemp = currentTime + minPosRoot(coef, Val(2))
            if timetemp < nextTime[i] 
                nextTime[i]=timetemp
            end
        else
            nextTime[i]=currentTime+1e-19 #
           println("**********************************nexttime[$i]==$currentTime")
            #= if  nextTime[1]==Inf &&  nextTime[2]==Inf
                println("equilibrium reached sys of 2 vars")
                @show currentTime
            end =#
        end

   #=coef=setindex(coef,q - x,1)
    time3 = currentTime + minPosRoot(coef, Val(2))
    time1 = time1 < time3 ? time1 : time3  =#  
   # nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    #= if q*q1<0 && abs(q) > 10.0*quantum[i] # uncomment did nothing
        time3=-(q/q1)-2*abs(quantum[i]/q1)
        if time3>0
        nextTime[i] = currentTime+time3 < nextTime[i] ? currentTime+time3 : nextTime[i]
        end
    end  =#   
   #=  else
        nextTime[i]=Inf
    end =#


end
#shortcut
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
   
   
    #= coef=setindex(coef,-quantum[i],1)
    coef=@SVector [#= q - x + =# -quantum[i], q1-x1,-x2]#
    time3 = currentTime + minPosRoot(coef, Val(2)) 
   # time3=time1 < time3 ? time1 : time3 
    coef=setindex(coef,quantum[i],1)
    time4 = currentTime + minPosRoot(coef, Val(2))
    nextTime[i] = time4 < time3 ? time4 : time3  
    #nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i] =#


    if nextTime[i] < currentTime  
        println("q neither heading towards x nor getting away by 2quantums")
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
   
    coef=@SVector [q - x , q1-x1,q2-x2,-x3]# x and q might get away even further(change of sign) and we want that to be no more than another quan
    nextTime[i] = currentTime + minPosRoot(coef, Val(3))
   

    if q-x >0
        coef=setindex(coef, -quantum[i],1)
        timetemp = currentTime + minPosRoot(coef, Val(3))
        if timetemp < nextTime[i] 
            nextTime[i]=timetemp
        end
    elseif  q-x <0
        coef=setindex(coef, quantum[i],1)
        timetemp = currentTime + minPosRoot(coef, Val(3))
        if timetemp < nextTime[i] 
            nextTime[i]=timetemp
        end
    else
        nextTime[i]=currentTime
        println("normally q and x at this point should not be equal")
    end

    
end


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

#= function Liqss_reComputeNextTime(::Val{3}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
    q2=qv[i][2]
    x3=xv[i][3]
    if  x3!=0.0 && a[i][i] != 0.0 # this is at the end of the step so my x is 1 quan away and derivatives coincide
        #nextTime[i]=currentTime+abs((q2-x2)/(x3))#4.247 ms (934 allocations: 64.02 KiB)
       #nextTime[i]=currentTime+sqrt(abs((x1)/(x3)))#654.351 μs (934 allocations: 64.02 KiB)#shifts up a little
       nextTime[i]=currentTime+cbrt(abs((q-x)/(x3)))#634.635 μs (934 allocations: 64.02 KiB)#*******************correct*******
     # nextTime[i]=currentTime+cbrt(abs((quantum[i])/(x3)))#632.743 μs (934 allocations: 64.02 KiB)#shifts up
    else
        nextTime[i]=Inf#a=0 means i dont depend on myself and it does not make sense to compute next based on me (x3)
    end
    coef=@SVector [q - x + 2*quantum[i], q1-x1,q2-x2,-x3]# x and q might get away even further(change of sign) and we want that to be no more than another quan
    time1 = currentTime + minPosRoot(coef, Val(3))
    coef=setindex(coef,q - x - 2*quantum[i],1)
    time2 = currentTime + minPosRoot(coef, Val(3))
    time1 = time1 < time2 ? time1 : time2    
    nextTime[i] = time1 < nextTime[i] ? time1 : nextTime[i]
    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end =#

















#######################################################################################################################################################
function updateLinearApprox(::Val{1},i::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=q[i][0]-qaux[i][1]
    if diffQ != 0.0
        a[i][i]=(x[i][1]-olddx[i][1])/diffQ
    else
        a[i][i]=0.0
    end
   #=  if i==2 && abs(a[2][2]+1.0)>1e-2
        println("prob in a22 which supposed to be 1")
        @show simt
        @show a
        @show u
        @show q[i][0],qaux[i][1]
        @show x[i][1],olddx[i][1]
        end
        if i==1 && abs(a[1][1]+1.0)>1e-2
          println("prob in a11 which supposed to be -1")
          @show simt
          @show a
          @show u
          @show q[i][0],qaux[i][1]
          @show x[i][1],olddx[i][1]
          end =#






    u[i][i][1]=x[i][1]-a[i][i]*q[i][0]   #if a==0 u same as derx meaning that in updateQ if derx> we dont have to check if u>0 ....note1
    return nothing
end
function updateLinearApprox(::Val{2},i::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=q[i][0]-qaux[i][1]
   # println("aii before updateoher= ",a[i][i])
    if diffQ != 0.0
        if abs(a[i][i])>1e-6
        a[i][i]=(x[i][1]-olddx[i][1])/diffQ
        end
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

#nothing
   #=  if i==2 && abs(a[2][2]+1.0)>1e-2
        println("prob in a22 which supposed to be 1")
        @show simt
        @show a
        @show u
        @show q[i][0],qaux[i][1]
        @show x[i][1],olddx[i][1]
        end
        if i==1 && abs(a[1][1]+1.0)>1e-2
          println("prob in a11 which supposed to be -1")
          @show simt
          @show a
          @show u
          @show q[i][0],qaux[i][1]
          @show x[i][1],olddx[i][1]
          end =#
   # println("aii before updateoher= ",a[i][i])
    u[i][i][1]=x[i][1]-a[i][i]*q[i][0]
    u[i][i][2]=2*x[i][2]-a[i][i]*q[i][1]
    tu[i]=simt  # comment did nothing but it makes sense to keep it because more accurate since u is changed
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
    tu[i]=simt  # uncomment did nothing
   
    return nothing
end
