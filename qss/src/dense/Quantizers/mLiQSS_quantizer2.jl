

function nmisCycle_and_simulUpdate(cacherealPosi,cacherealPosj,aij,aji,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{2},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{2,Float64}},qaux::Vector{MVector{2,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
  

  cacheA[1]=0.0;exacteA(q,d,cacheA,index,index)
   aii=cacheA[1]
   cacheA[1]=0.0; exacteA(q,d,cacheA,j,j)
   ajj=cacheA[1]


   uii=dxaux[index][1]-aii*qaux[index][1]
   ui2=dxaux[index][2]-aii*qaux[index][2]

  #=  uii=x[index][1]-aii*q[index][0]
   ui2=x[index][2]-aii*q[index][1] =#

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
    if DEBUG2  @show ddxj end
  end
  iscycle=false


    qjplus=xjaux-sign(ddxj)*quanj
    hj=sqrt(2*quanj/abs(ddxj))#
    α1=1-hj*ajj
    if abs(α1)==0.0
      α1=1e-30
      @show α1
    end
    dqjplus=(aji*(qi+hj*qi1)+ajj*qjplus+uji+hj*uji2)/α1
 

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
        if DEBUG2 @show ddxi end
      end
      hi=sqrt(2*quani/abs(ddxi))
      βidir=dxi+hi*ddxi/2
      βjdir=dxj+hj*ddxj/2
      βidth=dxithrow+hi*ddxithrow/2
      αidir=xi1+hi*xi2/2
      βj=xj1+hj*xj2/2
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
        coef_signig=100
          if (abs(dxi)>(abs(xi1)*coef_signig) #= && (abs(ddxi)>(abs(xi2)*coef_signig)||abs(xi2)>(abs(ddxi)*coef_signig)) =#) || (abs(xi1)>(abs(dxi)*coef_signig) #= &&  (abs(ddxi)>(abs(xi2)*coef_signig)||abs(xi2)>(abs(ddxi)*coef_signig)) =#)
            iscycle=true
          end
        #=  if abs(βidir-αidir)>(abs(βidir+αidir)/2)  
            iscycle=true
          end =#
   end


if (ddxithrow>0 && dxithrow<0 && qi1>0) || (ddxithrow<0 && dxithrow>0 && qi1<0)
  xmin=xi-dxithrow*dxithrow/ddxithrow
  if abs(xmin-xi)/abs(qi-xi)>0.8 # xmin is close to q and a change in q  might flip its sign of dq
    iscycle=false #
  end
end


#= if (ddxi>0 && dxi<0 && qi1>0) || (ddxi<0 && dxi>0 && qi1<0)
  xmin=xi-dxi*dxi/ddxi
  if abs(xmin-xi)/abs(qi-xi)>0.8 # xmin is close to q and a change in q  might flip its sign of dq
  
    iscycle=false #
  end

end =#

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

 #=    if DEBUG && 0.00047303128867631024 <=simt<=0.00047303497111002686 && (index==2 || index==1)
      println("*****check cycle conditions*********")
      println("i var: $index")
      @show xi1,dxithrow,dxi 
      @show xi2,ddxithrow,ddxi
      @show aii,qi1,aij,dqjplus,uij2
     #=  @show qi,hi
      @show βidir,βidth , βidir,dirI
      println("j var: $j")
      @show dxj,xj1,hj
      @show ddxj,xj2
      @show x[j],q[j]
      @show βjdir,βj ,dqjplus,newDiff =#
    end
 =#

     if iscycle
      trackSimul[1]+=1 
   


        h = ft-simt
        #h_=BigFloat(h)
        #aii,aij,aji,ajj,xi,xjaux,uij,uij2,uji,uji2=BigFloat(aii),BigFloat(aij),BigFloat(aji),BigFloat(ajj),BigFloat(xi),BigFloat(xjaux),BigFloat(uij),BigFloat(uij2),BigFloat(uji),BigFloat(uji2)
        qi,qj,Δ1=simulQ(quani,quanj,simt,aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2)
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
         # h_=BigFloat(h)
          qi,qj,Δ1=simulQ(quani,quanj,simt,aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2)
        end
        maxIter=10000
        while (abs(qi - xi) > 2.0*quani || abs(qj - xjaux) > 2.0*quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * sqrt(quani / abs(qi - xi));
          h2 = h * sqrt(quanj / abs(qj - xjaux));
          #=  h1 = h * (0.99*1.8*quani / abs(qi - xi));
          h2 = h * (0.99*1.8*quanj / abs(qj - xjaux)); =#
          h=min(h1,h2)
          #h_=BigFloat(h)
          qi,qj,Δ1=simulQ(quani,quanj,simt,aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2)
        end

        #qi,qj,Δ1=Float64(qi),Float64(qj),Float64(Δ1)
        if maxIter<9000  #ie simul step failed
          println("simulstep  $maxIter")
        
        end
        if maxIter==0  #ie simul step failed
          println("simulstep failed maxiter")
          return false
        end
        if  h<1e-20  #ie simul step failed
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
       

#= 
        if DEBUG && 0.0004710504575258516 <=simt<=0.0004740759390735583 && (index==2 || index==1)
          println("-------------end of simul------------")
          @show simt,h
          println("i var: $index")
          @show x[index],q[index],quani
          @show aij*qj+aii*qi+uij
          @show aij*q[j][1]+aii*q[index][1]+uij2
          println("j var: $j")
          @show x[j],q[j],quanj
          @show aji*qi+ajj*qj+uji
          @show aji*q[index][1]+ajj*q[j][1]+uji2
        end =#


      end #end if iscycle
    
    
  return iscycle
end


@inline function simulQ(quani,quanj,simt,aii::Float64,aij::Float64,aji::Float64,ajj::Float64,h::Float64,xi::Float64,xjaux::Float64,uij::Float64,uij2::Float64,uji::Float64,uji2::Float64)
#  @inline function simulQ(aii::BigFloat,aij::BigFloat,aji::BigFloat,ajj::BigFloat,h::BigFloat,xi::BigFloat,xjaux::BigFloat,uij::BigFloat,uij2::BigFloat,uji::BigFloat,uji2::BigFloat)
 
  h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
  aiijj=aii+ajj
  aiijj_=aij*aji-aii*ajj
  #Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
  Δ1=1.0-h*(aiijj)-h_2*(aiijj_)
  if abs(Δ1)==0.0
    Δ1=1e-30
    @show Δ1
  end
  
  Δ1_2=Δ1*Δ1
  Δ1_2_=1.0-2.0*h*aiijj+h_2*(aiijj*aiijj-2.0*aiijj_)+2.0*h_3*aiijj*aiijj_+h_4*aiijj_*aiijj_
  #= if abs(Δ1_2-Δ1_2_) >1e-3
    @show Δ1_2,Δ1_2_
  end =#
  αii=(aii*(1.0-h*ajj)+h*aij*aji)/Δ1
 # αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
  αij=aij/Δ1
  αji=aji/Δ1
  αjj=(h*aji*aij+(1.0-h*aii)*ajj)/Δ1
  βii_=1.0+h*(αii-aii)-h_2*(aii*αii+aij*αji)/2.0
  βii=1.0+(h_2*(aij*aji+aii*aii)/(2.0*Δ1)+h_3*aii*(aiijj_)/(2.0*Δ1))
  βij_=h*(αij-aij)-h_2*(aii*αij+aij*αjj)/2.0
  βij=h_2*aij*(aiijj)/(2.0*Δ1)+h_3*aij*( aiijj_)/(2.0*Δ1)
  βji_=h*(αji-aji)-h_2*(aji*αii+ajj*αji)/2.0
  βji=h_2*aji*(aiijj)/(2.0*Δ1)+h_3*aji*(aiijj_)/(2.0*Δ1)

  βjj_=1.0+h*(αjj-ajj)-h_2*(aji*αij+ajj*αjj)/2.0
  #β__jj=1+(h*h*(aji*aij+ajj*ajj)/2+h*h*h*(ajj*aij*aji-aii*ajj*ajj)/2)/Δ1
  βjj=1.0+(h_2*(aji*aij+ajj*ajj)/(2.0*Δ1)+h_3*ajj*( aiijj_)/(2.0*Δ1))
  

  Δ2__=βii*βjj-βij*βji
  
  
  Δ2__=1.0+(h_2*(aii*aii+ajj*ajj+2.0*aij*aji)+h_3*(aiijj)*(aiijj_))/(2.0*Δ1)+(h_4*(aiijj_)*(aiijj_)-h_5*(aiijj)*(aiijj_)*(aiijj_)-h_6*aiijj_*aiijj_*aiijj_)/(4.0*Δ1*Δ1)
  Δ2_=βii_*βjj_-βij_*βji_

  Δ2=(4.0-8.0*h*aiijj+h_2*(6.0*aiijj*aiijj-4.0*aiijj_)+h_3*(6.0*aiijj*aiijj_-2.0*aiijj*aiijj*aiijj)+h_4*(aiijj_*aiijj_-4.0*aiijj_*aiijj*aiijj)-3.0*h_5*aiijj*aiijj_*aiijj_-h_6*aiijj_*aiijj_*aiijj_)/(4.0*Δ1*Δ1) 
  if abs(Δ2-Δ2__) >10.0
    @show Δ2,Δ2__
    @show @show  aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2
  end
  if abs(Δ2)==0.0
    Δ2=1e-30
    @show Δ2
  end
  λii=(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij/2
  λij=(h*h*aii/2-h)*h*aij+h*h*aij*(1-h*aii)/2
  λji=h*h*aji/2*(1-h*ajj)+(h*h*ajj/2-h)*h*aji
  λjj=h*h*h*aij*aji/2+(h*h*ajj/2-h)*(1-h*aii)

  parti_=((λii*(uij+h*uij2)+λij*(uji+h*uji2))/Δ1)+(xi+h*uij+h*h*uij2/2)#part1[1]+xpart2[1]#
  parti=((xi-h*xi*(aiijj))-h_2*(aii*uij+aij*uji+uij2+2*xi*(aiijj_))/2+h_3*(-uij*(aiijj_)+uij2*ajj-aij*uji2)/2)/Δ1
  #parti=xi-h*h*((aii+ajj)*uij+aij*uji+uij2)/(2*Δ1)+h*h*h*(uij*(aii*ajj-aij*aji)+uij2*ajj-aij*uji2)/(2*Δ1)
 # parti=(-h*uij+h*h*(aii*uij-aij*uji+2*uij*ajj-2*uij2)/2+h*h*h*(aii*uij2-aii*ajj*uij+2*ajj*uij2+aji*aij*uij-aij*aii*uji-aij*uji2)/2+h*h*h*h*(aji*aij*uij2-aii*ajj*uij2)/2)/Δ1+(xi+h*uij+h*h*uij2/2)
# parti=((uij+h*uij2)*(h*h*aii/2-h)*(1-h*ajj)+h*h*h*aji*aij*(uij+h*uij2)/2+(h*h*aii/2-h)*h*aij*(uji+h*uji2)+h*h*aij*(1-h*aii)*(uji+h*uji2)/2)/Δ1+(xi+h*uij+h*h*uij2/2)
 #parti=(uij*h*h*aii/2-h*uij+h*h*h*aii*uij2/2-h*h*uij2-uij*h*h*h*aii*ajj/2+h*h*uij*ajj-h*h*h*h*ajj*aii*uij2/2+h*h*h*ajj*uij2+h*h*h*aji*aij*(uij+h*uij2)/2+(h*h*aii/2-h)*h*aij*(uji+h*uji2)+h*h*aij*(1-h*aii)*(uji+h*uji2)/2)/Δ1+(xi+h*uij+h*h*uij2/2)
 #parti_=(-h*uij+h*h*(-aij*uji/2+uij*aii/2-uij2+uij*ajj)+h*h*h*(-aij*uji2/2+aji*aij*uij/2+ajj*uij2+aii*uij2/2-uij*aii*ajj/2)+h*h*h*h*(aji*aij*uij2/2-ajj*aii*uij2/2))/Δ1+(h*uij+h*h*uij2/2)+xi#
 #parti=(-h*uij+h*h*(-aij*uji/2+uij*aii/2-uij2/2+uij*ajj)+h*h*h*(-aij*uji2/2+aji*aij*uij/2+ajj*uij2-ajj*uij2/2-uij*aii*ajj/2))/Δ1+h*uij+xi#
# if abs(parti-parti_)>1e-15
  #parti=parti_
    #@show parti,parti_,h,Δ1
    #@show  aii,aij,aji,ajj,h,xi,xjaux,uij,uij2,uji,uji2
   #=  @show (h*uij+h*h*uij2/2)
    @show uij,uij2 =#
  #end
  partj_=((λji*(uij+h*uij2)+λjj*(uji+h*uji2))/Δ1)+(xjaux+h*uji+h*h*uji2/2)#part1[2]+xpart2[2]#
  partj=((xjaux-h*xjaux*(aiijj))-h_2*((ajj)*uji+aji*uij+uji2+2*xjaux*(aiijj_))/2+h_3*(-uji*(aiijj_)+uji2*aii-aji*uij2)/2)/Δ1

  qi=((βjj/Δ2)*parti-(βij/Δ2)*partj)
  qj=((βii/Δ2)*partj-(βji/Δ2)*parti)
  qi_=((βjj_/Δ2_)*parti_-(βij_/Δ2_)*partj_)
  qj_=((βii_/Δ2_)*partj_-(βji_/Δ2_)*parti_)
 #=  if abs(qi-qi_)>1e-10
    @show qi,qi_,xi,quani,h
    @show simt
  end
  if abs(qj-qj_)>1e-10
    @show qj,qj_,xjaux,quanj,h
    @show simt
  end =#
 #=  if simt==11.163688670259043
    @show simt
    dif_=abs(xi-qi_);dif=abs(xi-qi)
    @show qi_,qi,xi,quani,dif_,dif
    difj_=abs(xjaux-qj_);difj=abs(xjaux-qj)
    @show qj_,qj,xjaux,quanj,difj_,difj
  end =#
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
      
     cacheA[1]=0.0; exacteA(qv,cacheA,i,i)
    
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