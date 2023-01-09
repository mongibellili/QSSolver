
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
function simulUpdate(::Val{1},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];uij=u[index][j][1];uji=u[j][index][1];xi=x[index][0];xj=x[j][0];x1i=x[index][1];x1j=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j]
  quani=quantum[index]
  qaux[j][1]=qj
  olddx[j][1]=x1j
  elapsed = simt - tx[j]
  xjaux = xj+elapsed*x1j
  dxj=aji*qi+ajj*qj+uji
  iscycle=false
  if dxj*x1j<0
    qjplus=xjaux+sign(dxj)*quanj
    dxi=aii*qi+aij*qjplus+uij
    if dxi*x1j<0
     # println("********simul update(val1) double if passed; simt= ",simt)
      iscycle=true  
      println("simul update happened at simt= ",simt)          
      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
      qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) 
        h1 = (abs(quani / x1i));h2 = (abs(quanj / x1j));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=100
      while (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (0.95*quani / abs(qi - xi));
        h2 = h * (0.95*quanj / abs(qj - xjaux));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      q[index][0]=qi# store back helper vars
      q[j][0]=qj
    end #end second dependecy check
 end # end outer dependency check
 return iscycle
end
         
function simulUpdate(::Val{2},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1]
  uii=u[index][index][1];uij=u[index][j][1];ujj=u[j][j][1];uji=u[j][index][1];uji2=u[j][index][2]
  quanj=quantum[j]
  quani=quantum[index]
  xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  e1 = simt - tx[j]
  xjaux = x[j](e1)# xAUX instead
   #e3=simt - tu[j]
   if debug
   println("----beginging of isCycle function: ")
   println("aii,aij,ajj,aji= a$index$index,a$index$j,a$j$j,a$j$index= ",aii," ",aij," ",ajj," ",aji)
   println("uii,uij,ujj,uji= u$index$index,u$index$j,u$j$j,u$j$index= ",uii," ",uij," ",ujj," ",uji)
   
   end
  
   # tx[j]=simt   # do not update tx[j] here ie do not uncomment because we are not doing a real update it is just a prediction
   # tu[j]=simt  # does not matter
   # tq[j]=simt   # does not matter
 
   e2 = simt - tq[j]
   qj=qj+e2*qj1  
  qaux[j][1]=qj
  #= q[j][0]=qj =#
  #@show xj1,e1,xi2,u[j][index][2]
  olddx[j][1]=xj1+e1*xj2#-e1*u[j][index][2]# # very huge spikes in many occasions
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
 # if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))    
    qjplus=xjaux-sign(ddxj)*quanj
    h=sqrt(quanj/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
    dqjplus=(aji*(qi+h*qi1)+ajj*qjplus+uji+h*uji2)/(1-h*ajj)
   ### 
    #u[index][j][1]=u[index][index][1]-a[index][j]*q[j][0]  # shifts down at 18
   
    u[index][j][1]=u[index][index][1]-aij*qaux[j][1]
    uij=u[index][j][1]
    u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
    uij2=u[index][j][2]
    dxi=aii*qi+aij*qjplus+uij
    ddxi=aii*qi1+aij*dqjplus+uij2
    if debug
    println("----in case the first if passed isCycle function modified these: ")
    @show uij
    @show uij2
    @show dxi,xi1
    @show ddxi,xi2
    end
  #  if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2))
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
        println("simul update happened at simt= ",simt) 
        A=[aii aij;aji ajj]
        I=[1 0;0 1]
        U=[uij;uji]
        U2=[uij2;uji2]
        X=[xi;xjaux]
        h1=0.0;h2=0.0
        h = ft-simt
        N=inv(I-h*A)
        Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2)        
        qi=Q[1]
        qj=Q[2]
        if (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) # removing this did nothing...check @btime later
          h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0
          h=min(h1,h2)
          N=inv(I-h*A)
          Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2)         
          qi=Q[1]
          qj=Q[2]
        end
        maxIter=600
        while (abs(qi - xi) > quani || abs(qj - xjaux) > quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * (0.95*quani / abs(qi - xi));
          h2 = h * (0.95*quanj / abs(qj - xjaux));
          h=min(h1,h2)
          N=inv(I-h*A)
          Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2)         
          qi=Q[1]
          qj=Q[2]
        end
       #=  if maxIter < 20
           @show maxIter
           @show simt
           @show a
        end =#
       if debug @show maxIter  end
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        Q1=N*(A*Q+U+h*U2)
        q[index][1]=Q1[1]# store back helper vars
        q[j][1]=Q1[2]

        if debug
          println("h used was = ",h) 
          @show q[index]
          @show q[j]       
          println("end of iscycle function")
         end
        tq[j]=simt
 #     end #end second dependecy check
#  end # end outer dependency check
  if debug  
    println("end of iscycle function")
   end
   nextTime[i]=simt+h
  return iscycle
end



#######################################################################################################################################################
#= function updateOtherApprox(::Val{1},j::Int,index::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
      diffQ=(q[index][0]-qaux[index][1])
    # @show q[index][0],qaux[index][1]
    # @show x[j][1],olddx[j][1]
    # println("aji before updateoher= ",a[j][index])
    #if 1.21<simt<1.226597
    #@show simt 
    #@show j
    # println("aji before updateoher= ",a[j][index])
    # @show q[index][0],qaux[index][1]
    # @show x[j][1],olddx[j][1]
      # end
        if diffQ!=0
        a[j][index]=(x[j][1]-olddx[j][1])/diffQ
        else
        a[j][index]=0.0
        end
      # if 1.21<simt<1.226597
      #   println("aji afterupdateoher= ",a[j][index])
          
        # end
      # @show a[j][j],a[j][index]
      # println("u inside updateOther before update= ",u[j][index])
    #    u[j][index][1]=x[j][1]-a[j][j]*q[j][0]-a[j][index]*q[index][0]
      # println("u inside updateOther after update= ",u[j][index])
        return nothing
end
function updateOtherApprox(::Val{2},j::Int,index::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
  diffQ=q[index][0]-qaux[index][1]
  # println("aji before updateoher= ",a[j][index])
  if diffQ != 0.0
      a[j][index]=(x[j][1]-olddx[j][1])/diffQ
  else
      a[j][index]=0.0
  end
  #  @show a[j][j],a[j][index]
  # if 1.54<simt<1.546597
  # println("a$j$i afterupdateoher= ",a[j][index])
    
  # end
  # println("u inside updateOther before update= ",u[j][index])
     u[j][index][1]=x[j][1]-a[j][j]*q[j][0]-a[j][index]*q[index][0]
   u[j][index][2]=2*x[j][2]-a[j][j]*q[j][1]-a[j][index]*q[index][1]
  #  println("u inside updateOther after update= ",u[j][index])
    #tu[index]=simt  # comment did nothing but it makes sense to keep it because more accurate since u is changed
  return nothing
end =#
#= function updateOlddx(::Val{1}, k::Int,x::Vector{Taylor0{Float64}},olddx::MVector{T,MVector{O,Float64}},elapsed::Float64)where{T,O}
  olddx[k][1]=x[k][1]
end
function updateOlddx(::Val{2}, k::Int,x::Vector{Taylor0{Float64}},olddx::MVector{T,MVector{O,Float64}},elapsed::Float64)where{T,O}
  olddx[k][1]=x[k][1]+2*x[k][2]*elapsed
end =#