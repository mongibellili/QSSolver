
#= function mupdateQ(::Val{1},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
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
      h = ft-simt
      q = (x + h * u) /(1 - h * a)
      if (abs(q - x) > 2 * quantum[i]) # removing this did nothing...check @btime later
        h = (abs(2*quantum[i] / dx));
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
end =#

#= function mupdateQ(::Val{2},i::Int, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},av::MVector{T,MVector{T,Float64}},uv::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
  q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; u1=uv[i][i][1]; u2=uv[i][i][2]
 #=  println("inside mupdate: before qaux[i][1]= ",qaux[i][1])
  qaux[i][1]=q#+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
  println("inside mupdate:after qaux[i][1]= ",qaux[i][1]) =#
  # println("qaux i inside updateQ= ",qaux[i][1])
  #q=qaux[i][1]# not needed...q used only in approx ddx which not needed to be exacte
 # q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; u1=uv[i][i][1]; u2=uv[i][i][2]

  #= if debug
    println("qaux$i was = ",qaux[i][1])
    println("q$i was = ",qv[i][0])
  end =#
  e=simt-tq[i]
  qaux[i][1]=q+e*q1
  if debug println("inside mupdate qaux$i becomes = ",qaux[i][1]) end
  qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent
  #tq[i]=simt
  #j=3-i    #this is for debugging only
  olddx[i][1]=x1#-e*uv[i][j][2]# 
  if debug println("olddx[i][1] becomes ",olddx[i][1]) end
  u1=u1+(simt-tu[i])*u2 # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
  uv[i][i][1]=u1
  tu[i]=simt  
 if debug 
  println("u$i$i= ",u1)
  println("tu$i= ",tu[i])
 end
 
  # olddx[i][2]=2*x2# 
  #=  if 14.1 <simt <=16.653630365539454
        @show simt,av
       
 end =#
  ddx=x2
  a=av[i][i]
  quan=quantum[i]
  h=0.0
  if a!=0.0
      if ddx ==0.0
          ddx=a*a*q+a*u1 +u2
          if ddx==0.0
              ddx=1e-26# changing -40 to -6 nothing changed
              println("ddx is zero")
          end
      end
      h = ft-simt
      q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
               (1 - h * a + h * h * a * a / 2)
      if (abs(q - x) > 2 * quan) # removing this did nothing...check @btime later
        h = sqrt(abs(2*quan / ddx)) # sqrt highly recommended...removing it leads to many sim steps..//2* is necessary in 2*quan when using ddx
        q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                 (1 - h * a + h * h * a * a / 2)
                # @show h
      end
      maxIter=4000
      while (abs(q - x) > 2 * quan) && (maxIter>0)
        maxIter-=1
        h = h *sqrt(2*quan / abs(q - x))
        q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                 (1 - h * a + h * h * a * a / 2)
      end
     #=  if maxIter < 200
       # println("maxiter of mpudate= ",maxIter)
        

      end =#
      q1=(a*q+u1+h*u2)/(1-h*a)  #later investigate 1=h*a
  else
      #ddx=u2
      if x2>0.0  #if ddx>0.0 same results
          q=x-quan
      else# elseif x2<0   ......else q=x??
          q=x+quan
      end
      if x2!=0.0
          h=sqrt(abs(2*quan/x2))   #sqrt necessary with u2
          q1=x1+h*x2  #(250 allocations: 16.07 KiB)
      else
          q1=x1
      end 
  end
  qv[i][0]=q
  qv[i][1]=q1  
 if debug 
  println("h used was = ",h) 
  println("inside mupdate q$i becomes = ",qv[i][0]) 
  println("inside mupdate derq$i becomes = ",qv[i][1]) 
  end
  return nothing
end =#

#############################################################################################################################
function isCycle_and_simulUpdate(::Val{1},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];uij=u[index][j][1];uji=u[j][index][1];xi=x[index][0];xj=x[j][0];x1i=x[index][1];x1j=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j]
  quani=quantum[index]
  qaux[j][1]=qj
  olddx[j][1]=x1j
 
  elapsed = simt - tx[j]
  x[j][0]= xj+elapsed*x1j
  xjaux=x[j][0]
  tx[j]=simt
  u[j][index][1]=u[j][j][1]-aji*qaux[index][1]
  uji=u[j][index][1]
  dxj=aji*qi+ajj*qj+uji
  iscycle=false
  h=0.0;h1=0.0;h2=0.0
  
  if dxj*x1j<0
    qjplus=xjaux+sign(dxj)*quanj
    u[index][j][1]=u[index][index][1]-aij*qj#qaux[j][1]
    uij=u[index][j][1]
    dxi=aii*qi+aij*qjplus+uij
    if dxi*x1i<0
     #=  qaux[j][1]=qj
      olddx[j][1]=x1j =#
     # println("********simul update(val1) double if passed; simt= ",simt)
      iscycle=true  
      println("***************************simul val1 update happened at simt= ",simt)   
      @show index       
      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
      qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
        h1 = (abs(quani / x1i));h2 = (abs(quanj / x1j));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (0.98*quani / abs(qi - xi));
        h2 = h * (0.98*quanj / abs(qj - xjaux));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
        if maxIter < 800
          println("maxiter of updateQ      = ",maxIter)
          @show  h, q-x
          
         end
      end

     #=  if 10.0>simt>7.676 && abs(x[2][0]-0.7)<3*quantum[2]
        @show h,h1,h2
        @show xi,qi,quani
        @show xj,qj,quanj
        @show a
          @show u
      end =#

      q[index][0]=qi# store back helper vars
      q[j][0]=qj
      tq[j]=simt

    end #end second dependecy check
 end # end outer dependency check
 return iscycle
end
         
function isCycle_and_simulUpdate(::Val{2},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1]
  uii=u[index][index][1];uij=u[index][j][1];
  ujj=u[j][j][1]#;uji=u[j][index][1]#;uji2=u[j][index][2]
  quanj=quantum[j]
  quani=quantum[index]
  xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  e1 = simt - tx[j]
  x[j][0]= x[j](e1)# xAUX instead
  xjaux=x[j][0]
  tx[j]=simt


  #= if abs(aii)>1e3
    @show index
  end =#
   #e3=simt - tu[j]
   #= if debug
   println("----beginging of isCycle function: ")
   println("aii,aij,ajj,aji= a$index$index,a$index$j,a$j$j,a$j$index= ",aii," ",aij," ",ajj," ",aji)
   println("uii,uij,ujj,uji= u$index$index,u$index$j,u$j$j,u$j$index= ",uii," ",uij," ",ujj," ",uji)
   
   end =#
  
   # tx[j]=simt   # do not update tx[j] here ie do not uncomment because we are not doing a real update it is just a prediction
   # tu[j]=simt  # does not matter
   # tq[j]=simt   # does not matter
 
   e2 = simt - tq[j]
   qj=qj+e2*qj1  
  qaux[j][1]=qj
  #= q[j][0]=qj =#
  #@show xj1,e1,xi2,u[j][index][2]
  xj1=x[j][1]+e1*xj2#-e1*u[j][index][2]# # very huge spikes in many occasions
  #olddx[j][1]=xj1
 # x[j][1]=xj1
  e3=simt - tu[j]
  ujj=ujj+e3*u[j][j][2]  #e3 does not work here
  tu[j]=simt

  u[j][index][1]=ujj-aji*qaux[index][1]# using q[i][0] creates a really huge bump at 18 (no go) because we want to elaps-update uji
  uji=u[j][index][1]
  
     u[j][index][2]=u[j][j][2]-a[j][index]*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
  # println("-------------------------a bug is here-----------------------------------------------")
 #u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # more cycles ...shaky with no bumps
  uji2=u[j][index][2]
  
  #@show uji2
  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2
  iscycle=false
  u[j][j][1]=ujj
  if debug
  println("----isCycle function modified these: ")
  @show xjaux
 
  @show qaux[j][1],olddx[j][1]
  
  @show  ujj,uji2
  
  @show dxj,xj1
  @show ddxj,xj2
  end
  #u[index][index][1]=u[index][index][1]+e1*u[index][index][2] # already updated in mupdate
  #= u[index][j][1]=u[index][index][1]-aij*qaux[j][1]
  uij=u[index][j][1]
  u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
    uij2=u[index][j][2] =#


    #= if abs(u[1][2][1])>1e-2 || abs(u[2][1][1]-2020)>1e-2 #sys2 testing
      @show simt
      @show index
      @show u[1][1][1],x[1][1]
      @show u[1][2][1]
      @show u[2][1][1]
      @show u[2][2][1],x[2][1]
    end =#



    u[index][j][1]=u[index][index][1]-aij*qaux[j][1]
    uij=u[index][j][1]
    u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
    uij2=u[index][j][2]


  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))    
    qjplus=xjaux-sign(ddxj)*quanj
    h=sqrt(quanj/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
    dqjplus=(aji*(qi+h*qi1)+ajj*qjplus+uji+h*uji2)/(1-h*ajj)
   ### 
    #u[index][j][1]=u[index][index][1]-a[index][j]*q[j][0]  # shifts down at 18
   
   #=  u[index][j][1]=u[index][index][1]-aij*qaux[j][1]
    uij=u[index][j][1]
    u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
    uij2=u[index][j][2] =#
    dxi=aii*qi+aij*qjplus+uij
    ddxi=aii*qi1+aij*dqjplus+uij2
    if debug
    println("----in case the first if passed isCycle function modified these: ")
    @show uij
    @show uij2
    @show dxi,xi1
    @show ddxi,xi2
    end
    if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2))
      #if abs(ddxi)>1000||abs(ddxj)>1000
       #=  @show simt
       @show dxi,xi1
       @show ddxi,xi2
       @show dxj,xj1
       @show ddxj,xj2
       @show a
       @show u =#
     # end
   #  println("********simul update(val2) double if passed; simt= ",simt)
        iscycle=true
       #= if simt==847.5459131958594
        println("*************simul val2 update happened at simt= ",simt) 
        @show index
       end =#
        A=[aii aij;aji ajj]
        I=[1.0 0.0;0.0 1.0]
        U=[uij;uji]
        U2=[uij2;uji2]
        X=[xi;xjaux]

        h = ft-simt
        N=inv(I-h*A)
        Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))    
       # Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2/2)      
        qi=Q[1]
        qj=Q[2]
        if (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) # removing this did nothing...check @btime later
          h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0
          h=min(h1,h2)
          N=inv(I-h*A)
          Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+(X+(h*U+(h*h*U2/2))))
         # Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2/2)  
          qi=Q[1]
          qj=Q[2]
        end
        maxIter=600
        while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * (1.96*quani / abs(qi - xi));
          h2 = h * (1.96*quanj / abs(qj - xjaux));
          h=min(h1,h2)
          N=inv(I-h*A)
         
          Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+( X+(h*U+(h*h*U2/2))) )
           
          
          qi=Q[1]
          qj=Q[2]
          
        end
        L=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+   X+ h*U+ h*h*U2/2    )
        #= if (abs(Q[1]-L[1])>1e-15 || abs(Q[2]-L[2])>1e-15) && 350<simt<350.999
          @show Q[1],X[1]
          @show Q[2],X[2]
          @show L
          @show A
          @show h
          @show simt
          @show U
          @show U2
          
         end =#
        #=  if 350.95<simt<350.96
          @show Q[1],X[1]
          @show Q[2],X[2]
          @show L
          @show A
          @show h
          @show simt
          @show U
          @show U2
          println("------------")
          
         end =#
        #= if simt>350
          @show simt
         end =#
        if maxIter < 20
           @show maxIter
           @show simt
           @show a
        end
       if debug @show maxIter  end
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        Q1=N*(A*Q+U+h*U2)
        q[index][1]=Q1[1]# store back helper vars
        q[j][1]=Q1[2]

        if  abs(qi+h*Q1[1]-xi)<1e-8  || abs(qj+h*Q1[2]-xjaux)<1e-8
          println("next x is about the same this x")
          @show h
          @show qi+h*Q1[1],xi
          @show qj+h*Q1[2],xjaux
          @show simt
          @show j
        end
        if abs(qi-xi)<1e-12 || abs(qj-xjaux)<1e-12
          println("could not throw q")
          @show qi,xi
          @show qj,xjaux
          @show simt
          @show j
        end
        if debug
          println("h used was = ",h) 
          @show q[index]
          @show q[j]       
          println("end of iscycle function")
         end
        tq[j]=simt
      end #end second dependecy check
  end # end outer dependency check
  if debug  
    println("end of iscycle function")
   end
  return iscycle
end


#######################################################################################################################################################
function Liqss_simulreComputeNextTime(::Val{1}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
  dt=0.0
  q=qv[i][0]
  x=xv[i][0]
  if xv[i][1] !=0.0  #&& abs(q-x)>quantum[i]/10
      if abs(q-x)>quantum[i]/10 # if q and x are very close then equilibrium no need to enter to compute next
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
                  println("liqss1 simulrecompute: x heading neither towards nor 2delta away! at time= ",currentTime) 
                  nextTime[i]=currentTime+1e-1
              end =#
          else#q=x equilibruim point reached#not needed with if abs(q-x)>quantum[i]/10
              nextTime[i]=Inf
              println("simulrecompval1 dt==0")
          end
          #= if 10.0>currentTime>7.676
              println(" nextTime[$i]= ", nextTime[i])
            end =#
      else
          nextTime[i]=Inf
          if nextTime[1]==Inf && nextTime[2]==Inf
          println("!!!!!!!!!!!!!!!!!!!!!!!!liqss1 simulrecompute: abs(q-x)<quantum[$i]/10!   equilibrium at time? ",currentTime) 
          end
      end
  else
      nextTime[i]=Inf
  end
end

#= function Liqss_simulreComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
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


    if nextTime[i]==Inf
        println("q neither heading towards x nor getting away by 2quantums")
       # nextTime[i]=currentTime+1e-6
    end



    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end =#

#= function Liqss_simulreComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
    q=qv[i][0]
    x=xv[i][0]
    q1=qv[i][1]
    x1=xv[i][1]
    x2=xv[i][2]
   #=  coef=@SVector [q-x, q1-x1,-x2]#
        nextTime[i]=currentTime + minPosRoot(coef, Val(2)) =#
        if q-x >0
           # coef=setindex(coef, -quantum[i],1)
            coef=@SVector [-quantum[i], q1-x1,-x2]
            timetemp = currentTime + minPosRoot(coef, Val(2))
            if timetemp < nextTime[i] 
                nextTime[i]=timetemp
            end
        elseif  q-x <0
            #coef=setindex(coef, quantum[i],1)
            coef=@SVector [quantum[i], q1-x1,-x2]
            timetemp = currentTime + minPosRoot(coef, Val(2))
            if timetemp < nextTime[i] 
                nextTime[i]=timetemp
            end
        else
            nextTime[i]=currentTime
            println("normally q and x at this point should not be equal")
        end



    #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
        time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
        nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
    end  =#   
end =#

function Liqss_simulreComputeNextTime(::Val{2}, i::Int, currentTime::Float64, nextTime::MVector{T,Float64}, xv::Vector{Taylor0{Float64}},qv::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}})where{T}
  q=qv[i][0]
  x=xv[i][0]
  q1=qv[i][1]
  x1=xv[i][1]
  x2=xv[i][2]
  if abs(q-x)>quantum[i]/10
  coef=@SVector [q-x, q1-x1,-x2]#
      nextTime[i]=currentTime + minPosRoot(coef, Val(2))
      if q-x >0
          coef=setindex(coef, -quantum[i],1)
          timetemp = currentTime + minPosRoot(coef, Val(2))
          if timetemp < nextTime[i] 
              nextTime[i]=timetemp
          end
      elseif  q-x <0
          coef=setindex(coef, quantum[i],1)
          timetemp = currentTime + minPosRoot(coef, Val(2))
          if timetemp < nextTime[i] 
              nextTime[i]=timetemp
          end
      else
          nextTime[i]=Inf #
          println("equilibrium reached")
      end
    else
      nextTime[i]=Inf
      if nextTime[1]==Inf && nextTime[2]==Inf
      println("!!!!!!!!!!!!!!!!!!!!!!!!liqss2 simulrecompute: abs(q-x)<quantum[$i]/10!   equilibrium at time? ",currentTime) 
      end
  end
  #= if q*q1<0 && a[i][i] > 10.0*quantum[i] # uncomment did nothing
      time3=currentTime-q/a[i][i]-2*abs(quantum[i]/q1)
      nextTime[i] = time3 < nextTime[i] ? time3 : nextTime[i]
  end  =#   
end


#######################################################################################################################################################
function updateOtherApprox(::Val{1},j::Int,index::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
      diffQ=(q[index][0]-qaux[index][1])
    # @show q[index][0],qaux[index][1]
    # @show x[j][1],olddx[j][1]
    # println("aji before updateoher= ",a[j][index])
                   #=  if 1.21<simt<1.426597
                    @show simt 
                    #@show j
                    println("aji before updateoher= ",a[j][index]) 
                    # @show q[index][0],qaux[index][1]
                    # @show x[j][1],olddx[j][1]
                      end=#
        if diffQ!=0
        a[j][index]=(x[j][1]-olddx[j][1])/diffQ
        else
        a[j][index]=0.0
        end
                    #= if j==2 && abs(a[2][1]-1.0)>1e-2
                      println("prob in a21 which supposed to be 1")
                      @show simt
                      @show a
                      @show u
                      @show q[index][0],qaux[index][1]
                      @show x[j][1],olddx[j][1]
                      end
                      if j==1 && abs(a[1][2]+1.0)>1e-2
                        println("prob in a12 which supposed to be -1")
                        @show simt
                        @show a
                        @show u
                        @show q[index][0],qaux[index][1]
                        @show x[j][1],olddx[j][1]
                        end =#
      # @show a[j][j],a[j][index]
      # println("u inside updateOther before update= ",u[j][index])
    #    u[j][index][1]=x[j][1]-a[j][j]*q[j][0]-a[j][index]*q[index][0]
        #= if j==2 && abs(a[2][1]-1.0)>1e-2
          println("if i fix u, it becomes: ")
          @show u
         
          end
          if j==1 && abs(a[1][2]+1.0)>1e-2
            println("if i fix u, it becomes: ")
            @show u
           
            end =#
      # println("u inside updateOther after update= ",u[j][index])
        return nothing
end
function updateOtherApprox(::Val{2},k::Int,j::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
  diffQ=q[j][0]-qaux[j][1]
  # println("aji before updateoher= ",a[j][index])
  if diffQ != 0.0
      a[k][j]=(x[k][1]-olddx[k][1])/diffQ
  else
      a[k][j]=0.0
  end
  #  @show a[j][j],a[j][index]
  # if 1.54<simt<1.546597
  # println("a$j$i afterupdateoher= ",a[j][index])
    
  # end
  # println("u inside updateOther before update= ",u[j][index])
   #=   u[k][j][1]=x[k][1]-a[k][k]*q[k][0]-a[k][j]*q[j][0]
   u[k][j][2]=2*x[k][2]-a[k][k]*q[k][1]-a[k][j]*q[j][1] =#
  #  println("u inside updateOther after update= ",u[j][index])
    #tu[index]=simt  # comment did nothing but it makes sense to keep it because more accurate since u is changed
  return nothing
end
#= function updateOlddx(::Val{1}, k::Int,x::Vector{Taylor0{Float64}},olddx::MVector{T,MVector{O,Float64}},elapsed::Float64)where{T,O}
  olddx[k][1]=x[k][1]
end
function updateOlddx(::Val{2}, k::Int,x::Vector{Taylor0{Float64}},olddx::MVector{T,MVector{O,Float64}},elapsed::Float64)where{T,O}
  olddx[k][1]=x[k][1]+2*x[k][2]*elapsed
end =#