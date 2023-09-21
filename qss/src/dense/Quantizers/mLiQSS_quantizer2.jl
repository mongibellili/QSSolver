

function nmisCycle_and_simulUpdate(cacherealPosi,cacherealPosj,aij,aji,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{2},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{2,Float64}},qaux::Vector{MVector{2,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
  

   exacteA(q,cacheA,index,index)
   aii=cacheA[1]
   exacteA(q,cacheA,j,j)
   ajj=cacheA[1]


   uii=dxaux[index][1]-aii*qaux[index][1]

   ui2=dxaux[index][2]-aii*qaux[index][2]



  xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1];xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  #uii=u[index][index][1];ujj=u[j][j][1]#;uij=u[index][j][1];uji=u[j][index][1]#;uji2=u[j][index][2]
  quanj=quantum[j];quani=quantum[index];
  e1 = simt - tx[j];e2 = simt - tq[j];#= e3=simt - tu[j];tu[j]=simt;  =#
  prevXj=x[j][0]
  x[j][0]= x[j](e1);xjaux=x[j][0];tx[j]=simt

  qj=qj+e2*qj1  ;qaux[j][1]=qj;tq[j] = simt    ;q[j][0]=qj  

  xj1=x[j][1]+e1*xj2;
  newDiff=(xjaux-prevXj)
  #= dirj=direction[j]
  if newDiff*dirj <0.0
    dirj=-dirj 
  elseif newDiff==0 && dirj!=0.0
    dirj=0.0  
  elseif newDiff!=0 && dirj==0.0
    dirj=newDiff
  else
  end          
  direction[j]=dirj =#
  #ujj=ujj+e1*u[j][j][2]  
  ujj=xj1-ajj*qj
  #u[j][j][1]=ujj
  uji=ujj-aji*qaux[index][1]# 
  #uji=u[j][index][1]
  uj2=xj2-ajj*qj1###################################################-----------------------
  uji2=uj2-aji*qaux[index][2]#
 #u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
 # uji2=u[j][index][2] 
  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2
  if abs(ddxj)==0.0
    ddxj=1e-30
    @show ddxj
  end
  iscycle=false
    qjplus=xjaux-sign(ddxj)*quanj
    h=sqrt(2*quanj/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
    α1=1-h*ajj
    if abs(α1)==0.0
      α1=1e-30
      @show α1
    end
    dqjplus=(aji*(qi+h*qi1)+ajj*qjplus+uji+h*uji2)/α1
    uij=uii-aij*qaux[j][1]
    # uij=u[index][j][1]
     uij2=ui2-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
       #β=dxi+sqrt(abs(ddxi)*quani/2)
      #h2=sqrt(2*quani/abs(ddxi))
    
      #uij2=u[index][j][2]
      dxithrow=aii*qi+aij*qj+uij
      ddxithrow=aii*qi1+aij*qj1+uij2
      dxi=aii*qi+aij*qjplus+uij
      ddxi=aii*qi1+aij*dqjplus+uij2
      if abs(ddxi)==0.0
        ddxi=1e-30
        @show ddxi
      end
      hi=sqrt(2*quani/abs(ddxi))
      βidir=dxi+hi*ddxi/2
      βjdir=dxj+h*ddxj/2
      βidth=dxithrow+hi*ddxithrow/2
      βj=xj1+h*xj2/2
########condition:Union 
  #=   if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || dqjplus*newDiff<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || βidir*dirI<0.0
        iscycle=true
      end
    end =#
 ########condition:Union i
  #=   if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || dqjplus*newDiff<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) || βidir*dirI<0.0
        iscycle=true
      end
    end =#
 ########condition:combineDer Union i
 if abs(βjdir-βj)>(abs(βjdir+βj)/2)  || dqjplus*newDiff<0.0 
    
  if abs(βidir-βidth)>(abs(βidir+βidth)/2)  || βidir*dirI<0.0
    iscycle=true
  end
end

########condition:kinda signif alone i
  #=   if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2)) # || dqjplus*newDiff<0.0 
    
      if (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) #|| βidir*dirI<0.0
        iscycle=true
      end
    end =#

     ########condition:cond1 
  #=    if  dqjplus*newDiff<0.0 #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if βidir*dirI<0.0
        iscycle=true
      end
    end

    if  dqjplus*newDiff<0.0 #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) 
        iscycle=true
      end
    end

    if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))   #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if  βidir*dirI<0.0
        iscycle=true
      end
    end =#
 


     ########condition:cond1 i
 #=     if  dqjplus*newDiff<0.0 #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if βidir*dirI<0.0
        iscycle=true
      end
    end

    if  dqjplus*newDiff<0.0 #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if (abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) || abs(ddxi-ddxithrow)>(abs(ddxi+ddxithrow)/2)) 
        iscycle=true
      end
    end

    if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))   #= || (dqjplus==0.0 && newDiff!=0.0)  =##(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
    
      if  βidir*dirI<0.0
        iscycle=true
      end
    end =#

   #=  if index==1 && j==4 && 0.1<simt<1.1
      iscycle=true
    end =#
     if iscycle
        h = ft-simt
      
        qi,qj,Δ1=simulQ(aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2)
        if (abs(qi - xi) > 2.0*quani || abs(qj - xjaux) > 2.0*quanj) 
          h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0
         # h1 = sqrt(abs(2*quani/ddxi));h2 = sqrt(abs(2*quanj/ddxj)); 
        #=  if xi2==0.0
          xi2=1e-30
          @show xi2
         end
         if xj2==0.0
          xj2=1e-30
          @show xj2
         end
         h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));  =#  #later add derderX =1e-12 when x2==0
          h=min(h1,h2)
          qi,qj,Δ1=simulQ(aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2)
        end
        maxIter=10000
        while (abs(qi - xi) > 2.0*quani || abs(qj - xjaux) > 2.0*quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * sqrt(quani / abs(qi - xi));
          h2 = h * sqrt(quanj / abs(qj - xjaux));
         #=   h1 = h * (0.99*1.8*quani / abs(qi - xi));
          h2 = h * (0.99*1.8*quanj / abs(qj - xjaux)); =#
          h=min(h1,h2)
          qi,qj,Δ1=simulQ(aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2)
        end
        if maxIter==0  #ie simul step failed
          println("simulstep failed maxiter")
          return false
        end
        if  h<1e-30  #ie simul step failed
          println("simulstep failed small h=",h)
          @show simt,maxIter
          return false
        end
        q[index][0]=qi# store back helper vars
        q[j][0]=qj     
        q1parti=aii*qi+aij*qj+uij+h*uij2
        q1partj=aji*qi+ajj*qj+uji+h*uji2
        if abs(q1parti)==0.0
          q1parti=1e-30*sign(q1parti)
          @show q1parti,simt
        end
        if abs(q1partj)==0.0
          q1partj=1e-30
          @show q1partj,simt
        end
        q[index][1]=((1-h*ajj)/Δ1)*q1parti+(h*aij/Δ1)*q1partj# store back helper vars
        q[j][1]=(h*aji/Δ1)*q1parti+((1-h*aii)/Δ1)*q1partj
      end #end if iscycle
   
  return iscycle
end


@inline function simulQ(aii::Float64,aij::Float64,aji::Float64,ajj::Float64,h::Float64,xi::Float64,xjaux::Float64,uij::Float64,uij2::Float64,uji::Float64,uji2::Float64)
  #use h_2=h*h/2
  Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
  if abs(Δ1)==0.0
    Δ1=1e-30
    @show Δ1
  end
  αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
  αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
  αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
  αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
  βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
  βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
  βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
  βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2

  Δ2=βii*βjj-βij*βji
  if abs(Δ2)==0.0
    Δ2=1e-30
    @show Δ2
  end
  λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
  λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
  λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
  λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)

  parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
  partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#

  qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
  qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
  return (qi,qj,Δ1)
end






####################################nliqss################################################################
#= 
function nisCycle_and_simulUpdate(::Val{2},index::Int,j::Int#= ,prevStepVal::Float64 =#,direction::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},map::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{O,Float64}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},olddxSpec::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,qminus::Vector{Float64})where{O}
  # @timeit "inside nmisCycle block1" begin
   #aii=getA(Val(Sparsity),cacheA,a,index,index,map);ajj=getA(Val(Sparsity),cacheA,a,j,j,map);aij=getA(Val(Sparsity),cacheA,a,index,j,map);aji=getA(Val(Sparsity),cacheA,a,j,index,map)
   
   #aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index]

   map(q,cacheA,index,index)
   aii=cacheA[1]
   map(q,cacheA,j,j)
    ajj=cacheA[1]
    map(q,cacheA,index,j)
    aij=cacheA[1]
    map(q,cacheA,j,index)
    aji=cacheA[1]

  xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1];xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  # uii=u[index][index][1];#ujj=u[j][j][1]#;uij=u[index][j][1];uji=u[j][index][1]#;uji2=u[j][index][2]
   uii=dxaux[index][1]-aii*qaux[index][1]
    # u[index][index][1]=uii
    ui2=dxaux[index][2]-aii*qaux[index][2]
   # u[index][index][2]=ui2
   # ui2=u[index][index][2]
   quanj=quantum[j];quani=quantum[index];
     
   e1 = simt - tx[j];e2 = simt - tq[j];#= e3=simt - tu[j];tu[j]=simt;  =#
   x[j][0]= x[j](e1);xjaux=x[j][0];tx[j]=simt
   qminus[j]=qj
   qj=qj+e2*qj1  ;qaux[j][1]=qj;tq[j] = simt    ;q[j][0]=qj  
 
  
 
   xj1=x[j][1]+e1*xj2;olddxSpec[j][1]=xj1;olddx[j][1]=xj1
 
  #=  newDiff=(xjaux-prevStepVal)
   dirj=direction[j]
   if newDiff*dirj <0.0
     dirj=-dirj
  
   elseif newDiff==0 && dirj!=0.0
     dirj=0.0
    
   elseif newDiff!=0 && dirj==0.0
     dirj=newDiff
   else
  
   end          
   direction[j]=dirj =#
 
   #ujj=ujj+e1*u[j][j][2]  
   ujj=xj1-ajj*qj
   #u[j][j][1]=ujj
   uji=ujj-aji*qaux[index][1]# 
  # u[j][index][1]=uji
 
   #u[j][j][2]=xj2-ajj*qj1###################################################-----------------------
   uj2=xj2-ajj*qj1
  # u[j][j][2]=uj2
   uji2=uj2-aji*qaux[index][2]#
  #u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
  # u[j][index][2]=uji2
 
 
   #@show uji2
   dxj=aji*qi+ajj*qaux[j][1]+uji
   ddxj=aji*qi1+ajj*qj1+uji2
   #@show aji,ajj,uji2
   iscycle=false
   
    
 
     qjplus=xjaux-sign(ddxj)*quanj
     h=sqrt(2*quanj/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
     dqjplus=(aji*(qi+h*qi1)+ajj*qjplus+uji+h*uji2)/(1-h*ajj)
    
     
   #end#end init iscycle
 
     if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  #|| (dqjplus)*dirj<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
       #β=dxi+sqrt(abs(ddxi)*quani/2)
       #h2=sqrt(2*quani/abs(ddxi))
       uij=uii-aij*qaux[j][1]
     #  u[index][j][1]=uij
       uij2=ui2-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
      # u[index][j][2]=uij2
 
       dxi=aii*qi+aij*qjplus+uij
     ddxi=aii*qi1+aij*dqjplus+uij2
 
       βidir=dxi+sqrt(2*quani/abs(ddxi))*ddxi/2
      
       if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) #|| βidir*direction[index]<0.0
 
     
       #  @timeit "inside nmisCycle block2" begin
         iscycle=true
         h = ft-simt
 
         Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
         αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
        αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
        αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
        αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
        βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
        βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
        βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
        βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
 
        Δ2=βii*βjj-βij*βji
 
         λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
         λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
         λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
         λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)
 
         parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
         partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
 
         qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
          qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
 
         if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
           h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0
           h=min(h1,h2)
 
           Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
 
         αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
         αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
         αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
         αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
 
         βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
         βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
         βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
         βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
 
          Δ2=βii*βjj-βij*βji
  
         λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
         λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
         λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
         λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)
 
         parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
         partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
   
         
          qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
          qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
   
         end
        
         maxIter=600
         while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
           maxIter-=1
          #=  h1 = h * (0.98*quani / abs(qi - xi));
           h2 = h * (0.98*quanj / abs(qj - xjaux)); =#
           h1 = h * sqrt(quani / abs(qi - xi));
          h2 = h * sqrt(quanj / abs(qj - xjaux));
           h=min(h1,h2)
 
           Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
 
           αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
           αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
           αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
           αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
 
           βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
           βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
           βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
           βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
  
            Δ2=βii*βjj-βij*βji
   
 
 
         λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
         λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
         λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
         λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)
 
         parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
         partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
       
         qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
          qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
 
 
        
           
         end
     
         if maxIter < 1
            @show maxIter
            @show simt
            @show a
         end
 
 
      
         q[index][0]=qi# store back helper vars
         q[j][0]=qj
       
          q1parti=aii*qi+aij*qj+uij+h*uij2
          q1partj=aji*qi+ajj*qj+uji+h*uji2
       
          q[index][1]=((1-h*ajj)/Δ1)*q1parti+(h*aij/Δ1)*q1partj# store back helper vars
         q[j][1]=(h*aji/Δ1)*q1parti+((1-h*aii)/Δ1)*q1partj
      # end#end block2?????????????????????????????????????????????????????????????????
       end #end second dependecy check
     end # end outer dependency check
   return iscycle
end
 =#

#####################################old mliqss

function isCycle_and_simulUpdate(::Val{2},index::Int,j::Int,#= direction::Vector{Float64}, =# x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)where{O}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];
 # aii=getA(Val(Sparsity),cacheA,a,index,index,map);ajj=getA(Val(Sparsity),cacheA,a,j,j,map);aij=getA(Val(Sparsity),cacheA,a,index,j,map);aji=getA(Val(Sparsity),cacheA,a,j,index,map)
  xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1];xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  uii=u[index][index][1];ujj=u[j][j][1]#;uij=u[index][j][1];uji=u[j][index][1]#;uji2=u[j][index][2]
  quanj=quantum[j];quani=quantum[index]; 
  e1 = simt - tx[j];e2 = simt - tq[j];#= e3=simt - tu[j];tu[j]=simt;  =#
  x[j][0]= x[j](e1);xjaux=x[j][0];tx[j]=simt
   qj=qj+e2*qj1  ;qaux[j][1]=qj;tq[j] = simt    ;q[j][0]=qj  
  xj1=x[j][1]+e1*xj2#= ;olddxSpec[j][1]=xj1 =#;olddx[j][1]=xj1
 # ujj=ujj+e1*u[j][j][2]  
 ujj=xj1-ajj*qj
  u[j][j][1]=ujj
  u[j][index][1]=ujj-aji*qaux[index][1]# using q[i][0] creates a really huge bump at 18 (no go) because we want to elaps-update uji
  uji=u[j][index][1]
  u[j][j][2]=xj2-ajj*qj1###################################################-----------------------
  u[j][index][2]=u[j][j][2]-aji*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
 #u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
  uji2=u[j][index][2] 

  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2

  iscycle=false
  if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))
    qjplus=xjaux-sign(ddxj)*quanj
    h=sqrt(2*quanj/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
    dqjplus=(aji*(qi+h*qi1)+ajj*qjplus+uji+h*uji2)/(1-h*ajj)
    u[index][j][1]=uii-aij*qaux[j][1]
    uij=u[index][j][1]
    u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
    uij2=u[index][j][2]
    dxi=aii*qi+aij*qjplus+uij
    ddxi=aii*qi1+aij*dqjplus+uij2
  if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2))

        iscycle=true

            h = ft-simt

        Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
       αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
       αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
       αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1
       βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
       βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
       βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
       βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2

       Δ2=βii*βjj-βij*βji

        λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
        λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
        λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
        λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)

        parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
        partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#

        qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
         qj=((βii/Δ2)*partj-(βji/Δ2)*parti)

        if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
          h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0
          h=min(h1,h2)

          Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji

        αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
        αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
        αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
        αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1

        βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
        βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
        βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
        βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2

         Δ2=βii*βjj-βij*βji
 
        λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
        λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
        λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
        λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)

        parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
        partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
  
        
         qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
         qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
  
        end
       
        maxIter=600
        while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
          maxIter-=1
          #= h1 = h * (0.98*quani / abs(qi - xi));
          h2 = h * (0.98*quanj / abs(qj - xjaux)); =#
          h1 = h * sqrt(quani / abs(qi - xi));
          h2 = h * sqrt(quanj / abs(qj - xjaux));
          h=min(h1,h2)

          Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji

          αii=(aii*(1-h*ajj)+h*aij*aji)/Δ1
          αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
          αji=(aji*aii*h+(1-h*aii)*aji)/Δ1
          αjj=(h*aji*aij+(1-h*aii)*ajj)/Δ1

          βii=1+h*(αii-aii)-h*h*(aii*αii+aij*αji)/2
          βij=h*(αij-aij)-h*h*(aii*αij+aij*αjj)/2
          βji=h*(αji-aji)-h*h*(aji*αii+ajj*αji)/2
          βjj=1+h*(αjj-ajj)-h*h*(aji*αij+ajj*αjj)/2
 
           Δ2=βii*βjj-βij*βji
  


        λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
        λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
        λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
        λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)

        parti=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
        partj=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
      
        qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
         qj=((βii/Δ2)*partj-(βji/Δ2)*parti)


       
          
        end
    
        if maxIter < 1
           @show maxIter
           @show simt
           @show a
        end


     
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
      
         q1parti=aii*qi+aij*qj+uij+h*uij2
         q1partj=aji*qi+ajj*qj+uji+h*uji2
      
         q[index][1]=((1-h*ajj)/Δ1)*q1parti+(h*aij/Δ1)*q1partj# store back helper vars
        q[j][1]=(h*aji/Δ1)*q1parti+((1-h*aii)/Δ1)*q1partj
           
          
      
      end #end second dependecy check
    end # end outer dependency check
  return iscycle
end
#= 
function updateQ(::Val{2},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{O,Float64}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextTime::Vector{Float64})where{O}
      #a=getA(Val(Sparsity),cacheA,av,i,i,map)
      #a=av[i][i]
      
      exacteA(qv,cacheA,i,i)
    
        a=cacheA[1]
    # @show a,i
      q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; #u1=uv[i][i][1]; u2=uv[i][i][2]
      qaux[i][1]=q#+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
      qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent
      olddx[i][1]=x1#appears only here...updated here and used in updateApprox   
      olddx[i][2]=x2
      #u1=u1+(simt-tu[i])*u2 # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
      u1=x1-a*qaux[i][1]
      # uv[i][i][1]=u1
      dxaux[i][1]=x1
      dxaux[i][2]=x2
      u2=x2-a*q1
      # uv[i][i][2]= u2
      # tu[i]=simt  
        # olddx[i][2]=2*x2# 
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
            h = ft-simt
            #tempH1=h
            #q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /(1 - h * a + h * h * a * a / 2)
                    q=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
                  
            
            if (abs(q - x) > 2* quan) # removing this did nothing...check @btime later
              h = sqrt(abs(2*quan / ddx)) # sqrt highly recommended...removing it leads to many sim steps..//2* is necessary in 2*quan when using ddx
              q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                      (1 - h * a + h * h * a * a / 2)

              #  q = (x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
                    
            end
            maxIter=1000
          # tempH=h
        #=   if  (abs(q - x) >2*  quan)
          coef=@SVector [quan, a*quan,-(a*a*(x-quan)+a*u1+u2)/2]#
              h1= minPosRoot(coef, Val(2))
            coef=@SVector [-quan, -a*quan,-(a*a*(x+quan)+a*u1+u2)/2]#
              h2= minPosRoot(coef, Val(2))

            if h1<h2
                h=h1;q=x-quan
            else
                h=h2;q=x+quan
            end

        end =#

            while (abs(q - x) >2*  quan) && (maxIter>0) && (h>0)
                
              h = h *sqrt(quan / abs(q - x))
              q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                      (1 - h * a + h * h * a * a / 2)
                    #  q = (x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)       
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
      nextTime[i]=simt+h

        return nothing
 
end =#