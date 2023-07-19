
      

function nmisCycle_and_simulUpdate(::Val{2},prtyp::Val{Sparsity} ,cacheA::MVector{1,Int},map::Function ,index::Int,j::Int,prevStepVal::Float64,direction::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},olddxSpec::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,qminus::Vector{Float64})where{Sparsity,O}
  # @timeit "inside nmisCycle block1" begin
  aii=getA(Val(Sparsity),cacheA,a,index,index,map);ajj=getA(Val(Sparsity),cacheA,a,j,j,map);aij=getA(Val(Sparsity),cacheA,a,index,j,map);aji=getA(Val(Sparsity),cacheA,a,j,index,map)
  # aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index]
  xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1];xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  uii=u[index][index][1];ujj=u[j][j][1]#;uij=u[index][j][1];uji=u[j][index][1]#;uji2=u[j][index][2]
  quanj=quantum[j];quani=quantum[index];
  e1 = simt - tx[j];e2 = simt - tq[j];#= e3=simt - tu[j];tu[j]=simt;  =#
  x[j][0]= x[j](e1);xjaux=x[j][0];tx[j]=simt
  qminus[j]=qj
  qj=qj+e2*qj1  ;qaux[j][1]=qj;tq[j] = simt    ;q[j][0]=qj  

  xj1=x[j][1]+e1*xj2;olddxSpec[j][1]=xj1;olddx[j][1]=xj1
  newDiff=(xjaux-prevStepVal)
  dirj=direction[j]
  if newDiff*dirj <0.0
    dirj=-dirj 
  elseif newDiff==0 && dirj!=0.0
    dirj=0.0  
  elseif newDiff!=0 && dirj==0.0
    dirj=newDiff
  else
  end          
  direction[j]=dirj
  #ujj=ujj+e1*u[j][j][2]  
  ujj=xj1-ajj*qj
  u[j][j][1]=ujj
  u[j][index][1]=ujj-aji*qaux[index][1]# 
  uji=u[j][index][1]
  u[j][j][2]=xj2-ajj*qj1###################################################-----------------------
  u[j][index][2]=u[j][j][2]-aji*qaux[index][2]#
 #u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
  uji2=u[j][index][2] 
  dxj=aji*qi+ajj*qaux[j][1]+uji
  ddxj=aji*qi1+ajj*qj1+uji2
  iscycle=false
    qjplus=xjaux-sign(ddxj)*quanj
    h=sqrt(2*quanj/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
    dqjplus=(aji*(qi+h*qi1)+ajj*qjplus+uji+h*uji2)/(1-h*ajj)
    if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-xj2)>(abs(ddxj+xj2)/2))  || (dqjplus)*dirj<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
      #β=dxi+sqrt(abs(ddxi)*quani/2)
      #h2=sqrt(2*quani/abs(ddxi))
      u[index][j][1]=uii-aij*qaux[j][1]
      uij=u[index][j][1]
      u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
      uij2=u[index][j][2]
      dxi=aii*qi+aij*qjplus+uij
      ddxi=aii*qi1+aij*dqjplus+uij2
      βidir=dxi+sqrt(2*quani/abs(ddxi))*ddxi/2
      if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2)) || βidir*direction[index]<0.0
        iscycle=true
        h = ft-simt
        qi,qj,Δ1=simulQ(aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2)
        if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
          h1 = sqrt(abs(2*quani/xi2));h2 = sqrt(abs(2*quanj/xj2));   #later add derderX =1e-12 when x2==0
          h=min(h1,h2)
          qi,qj,Δ1=simulQ(aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2)
        end
        maxIter=600
        while (abs(qi - xi) >2* quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * sqrt(quani / abs(qi - xi));
          h2 = h * sqrt(quanj / abs(qj - xjaux));
          # h1 = h * (0.98*quani / abs(qi - xi));
         
        #  h2 = h * (0.98*quanj / abs(qj - xjaux));
          h=min(h1,h2)
          qi,qj,Δ1=simulQ(aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2)
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


@inline function simulQ(aii::Float64,aij::Float64,aji::Float64,ajj::Float64,h::Float64,xi::Float64,xjaux::Float64,uij::Float64,uij2::Float64,uji::Float64,uji2::Float64)
  #use h_2=h*h/2
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
  return (qi,qj,Δ1)
end

function updateOtherApprox(::Val{2},sparsity::Val{Sparsity},cacheA::MVector{1,Int},map::Function ,k::Int,j::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},simt::Float64)where{Sparsity,O}
  diffQ=q[j][0]-qaux[j][1]
  #akk=getA(Val(Sparsity),cacheA,a,k,k,map)
  #akk=a[k][k]
  if diffQ != 0.0
    #a[k][j]=(x[k][1]-olddx[k][1])/diffQ
    akjvalue=(x[k][1]-olddx[k][1])/diffQ
    #a[k][j]=akjvalue
    setA(Val(Sparsity),cacheA,a,k,j,map,akjvalue)
   # akjvalue=getA(Val(Sparsity),cacheA,a,k,j,map)
  #=   u[k][j][1]=x[k][1]-akk*q[k][0]-akjvalue*q[j][0]
    u[k][j][2]=2*x[k][2]-akk*q[k][1]-akjvalue*q[j][1] =#
  else
    setA(Val(Sparsity),cacheA,a,k,j,map,0.0)
   # a[k][j]=0.0
   #=  u[k][j][1]=x[k][1]-akk*q[k][0]
    u[k][j][2]=2*x[k][2]-akk*q[k][1] =#
  end
 
  
  return nothing
end
#= function updateOtherApprox2(::Val{2},sparsity::Val{Sparsity},cacheA::MVector{1,Int},map::Function,k::Int,j::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},simt::Float64)where{Sparsity,O}
  diffQ=q[j][0]-qaux[j][1]
  akk=getA(Val(Sparsity),cacheA,a,k,k,map)
  if diffQ != 0.0
    akjvalue=(x[k][1]-olddx[k][1])/diffQ
    setA(Val(Sparsity),cacheA,a,k,j,map,akjvalue)
   # akjvalue=getA(Val(Sparsity),cacheA,a,k,j,map)
    u[k][j][1]=x[k][1]-akk*q[k][0]-akjvalue*q[j][0]
    u[k][j][2]=2*x[k][2]-akk*q[k][1]-akjvalue*q[j][1]
  else
    setA(Val(Sparsity),cacheA,a,k,j,map,0.0)
    u[k][j][1]=x[k][1]-akk*q[k][0]
    u[k][j][2]=2*x[k][2]-akk*q[k][1]
  end
 
  
  return nothing
end =#
#= function updateOtherApprox(::Val{2},sparsity::Val{false},cacheA::MVector{1,Int},map::Function,k::Int,j::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},simt::Float64)where{T,O}
  diffQ=q[j][0]-qaux[j][1]
  #akk=getA(Val(Sparsity),cacheA,a,k,k,map)
  akk=a[k][k]
  if diffQ != 0.0
    akjvalue=(x[k][1]-olddx[k][1])/diffQ
  a[k][j]=akjvalue
  #=   if simt == 0.004395600232045285
      @show  akjvalue
      @show q[j][0]
      #@show olddx[k][1]
    end =#
   # setA(Val(Sparsity),cacheA,a,k,j,map,akjvalue)
   # akjvalue=getA(Val(Sparsity),cacheA,a,k,j,map)
    u[k][j][1]=x[k][1]-akk*q[k][0]-akjvalue*q[j][0]
    u[k][j][2]=2*x[k][2]-akk*q[k][1]-akjvalue*q[j][1]
  else
   # setA(Val(Sparsity),cacheA,a,k,j,map,0.0)
    a[k][j]=0.0
    u[k][j][1]=x[k][1]-akk*q[k][0]
    u[k][j][2]=2*x[k][2]-akk*q[k][1]
  end
 
  #=  u[k][j][1]=x[k][1]-a[k][k]*q[k][0]-a[k][j]*q[j][0]
   u[k][j][2]=2*x[k][2]-a[k][k]*q[k][1]-a[k][j]*q[j][1] =#
  #  println("u inside updateOther after update= ",u[j][index])
    #tu[index]=simt  # comment did nothing but it makes sense to keep it because more accurate since u is changed
  return nothing
end =#



####################################nliqss################################################################

function nisCycle_and_simulUpdate(::Val{2},prtyp::Val{Sparsity} ,cacheA::MVector{1,Int},map::Function ,index::Int,j::Int#= ,prevStepVal::Float64 =#,direction::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},olddxSpec::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,qminus::Vector{Float64})where{Sparsity,O}
  # @timeit "inside nmisCycle block1" begin
   aii=getA(Val(Sparsity),cacheA,a,index,index,map);ajj=getA(Val(Sparsity),cacheA,a,j,j,map);aij=getA(Val(Sparsity),cacheA,a,index,j,map);aji=getA(Val(Sparsity),cacheA,a,j,index,map)
   
   #aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index]
  xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1];xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
   uii=u[index][index][1];ujj=u[j][j][1]#;uij=u[index][j][1];uji=u[j][index][1]#;uji2=u[j][index][2]
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
   u[j][j][1]=ujj
   u[j][index][1]=ujj-aji*qaux[index][1]# 
   uji=u[j][index][1]
 
   u[j][j][2]=xj2-ajj*qj1###################################################-----------------------
 
   u[j][index][2]=u[j][j][2]-aji*qaux[index][2]#
  #u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # from article p20 line25 more cycles ...shaky with no bumps
   uji2=u[j][index][2] 
 
 
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
       u[index][j][1]=uii-aij*qaux[j][1]
       uij=u[index][j][1]
       u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
       uij2=u[index][j][2]
 
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


#####################################old mliqss

function isCycle_and_simulUpdate(::Val{2},prtyp::Val{Sparsity} ,cacheA::MVector{1,Int},map::Function ,index::Int,j::Int,#= direction::Vector{Float64}, =# x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},a::Vector{Vector{Float64}},u::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)where{Sparsity,O}
 #aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];
  aii=getA(Val(Sparsity),cacheA,a,index,index,map);ajj=getA(Val(Sparsity),cacheA,a,j,j,map);aij=getA(Val(Sparsity),cacheA,a,index,j,map);aji=getA(Val(Sparsity),cacheA,a,j,index,map)
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