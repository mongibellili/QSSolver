
#############################################################################################################################
#= function nmisCycle_and_simulUpdate(::Val{1},index::Int,j::Int,prevStepVal::Float64,direction::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},olddx::Vector{MVector{1,Float64}},olddxSpec::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,qminus::Vector{Float64})
  
  exacteA(q,cacheA,index,index)
  aii=cacheA[1]
  exacteA(q,cacheA,j,j)
  ajj=cacheA[1]
  exacteA(q,cacheA,index,j)
  aij=cacheA[1]
  exacteA(q,cacheA,j,index)
  aji=cacheA[1]

  uii=dxaux[index][1]-aii*qaux[index][1]
  #aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index]
 # uij=u[index][j][1];uji=u[j][index][1];
  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;olddx[j][1]=ẋj
 
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;xj=x[j][0]
  tx[j]=simt
  ujj=ẋj-ajj*qj
   
  uji=ujj-aji*qaux[index][1]
 
  dxj=aji*qi+ajj*qj+uji
  iscycle=false 
  if dxj*ẋj<0
    qjplus=xj+sign(dxj)*quanj
    uij=uii-aij*qj#qaux[j][1]
   
    dxi=aii*qi+aij*qjplus+uij
    if dxi*ẋi<0
      iscycle=true  
        
      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
      qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) 
        h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (1.96*quani / abs(qi - xi));
        h2 = h * (1.96*quanj / abs(qj - xj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          #println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        if maxIter < 1
          println("maxiter of updateQ      = ",maxIter)
         end
        end

      q[index][0]=qi# store back helper vars
      q[j][0]=qj
      tq[j]=simt

    end #end second dependecy check
 end # end outer dependency check
 return iscycle
end =#
 #this for old mliqss        


#= function nmisCycle_and_simulUpdate(::Val{1},index::Int,j::Int,prevStepVal::Float64,direction::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},olddx::Vector{MVector{1,Float64}},olddxSpec::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,qminus::Vector{Float64})
  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;olddx[j][1]=ẋj
 
  exacteA(q,cacheA,j,index)
  aji=cacheA[1]
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;
  xj=x[j][0]
  tx[j]=simt
  #dxj=aji*qi+ajj*qj+uji 
  iscycle=false 
  if (aji*(qi-qaux[index][1])+ẋj)*ẋj<0.0
    exacteA(q,cacheA,index,j)
    aij=cacheA[1]
    exacteA(q,cacheA,index,index)
    aii=cacheA[1]
   # qjplus=xj-sign(ẋj)*quanj
  #  dxi=aii*qi+aij*qjplus+uij
  # if dxi*ẋi<0.0
  # if ẋi*(ẋi+aii*(qi-qaux[index][1])+aij*(xj-qj-sign(ẋj)*quanj))<0.0
  if ẋi*(ẋi+aii*(qi-qaux[index][1])+aij*(-2*sign(ẋj)*quanj))<0.0
      iscycle=true  
      
      exacteA(q,cacheA,j,j)
      ajj=cacheA[1]
      
      ujj=ẋj-ajj*qj
      uji=ujj-aji*qaux[index][1]
      uii=olddx[index][1]-aii*qaux[index][1]
      uij=uii-aij*qj#qaux[j][1]
   


      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
      qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) 
        h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (1.96*quani / abs(qi - xi));
        h2 = h * (1.96*quanj / abs(qj - xj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          #println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        if maxIter < 1
          println("maxiter of updateQ      = ",maxIter)
         end
        end

      q[index][0]=qi# store back helper vars
      q[j][0]=qj
      tq[j]=simt

    end #end second dependecy check
 end # end outer dependency check
 return iscycle
end =#

#function similar to version that calculates a manually
#= function nmisCycle_and_simulUpdate(::Val{1},index::Int,j::Int,prevStepVal::Float64,direction::Vector{Float64}, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},olddx::Vector{MVector{1,Float64}},olddxSpec::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64,qminus::Vector{Float64})
  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;olddx[j][1]=ẋj
 
  exacteA(q,cacheA,j,index)
  aji=cacheA[1]
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;
  xj=x[j][0]
  tx[j]=simt
 
  exacteA(q,cacheA,index,j)
    aij=cacheA[1]
    exacteA(q,cacheA,index,index)
    aii=cacheA[1]
  exacteA(q,cacheA,j,j)
  ajj=cacheA[1]
  
  ujj=ẋj-ajj*qj
  uji=ujj-aji*qaux[index][1]
  uii=olddx[index][1]-aii*qaux[index][1]
  uij=uii-aij*qj#qaux[j][1]

  dxj=aji*qi+ajj*qaux[j][1]+uji
  iscycle=false
  #if (aji*(qi-qaux[index][1])+ẋj)*ẋj<0.0
  if (dxj*ẋj)<0.0 
   
    qjplus=xj+sign(dxj)*quanj
    dxi=aii*qi+aij*qjplus+uij
  if (dxi*ẋi)<0.0
  #if ẋi*(ẋi+aii*(qi-qaux[index][1])+aij*(-2*sign(ẋj)*quanj))<0.0
      iscycle=true  
      
   
   


      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
      qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) 
        h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (1.96*quani / abs(qi - xi));
        h2 = h * (1.96*quanj / abs(qj - xj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          #println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        if maxIter < 1
          println("maxiter of updateQ      = ",maxIter)
         end
        end

      q[index][0]=qi# store back helper vars
      q[j][0]=qj
      tq[j]=simt

    end #end second dependecy check
 end # end outer dependency check
 return iscycle
end =#

#= function nmisCycle_and_simulUpdate(::Val{1},index::Int,j::Int,dirI::Float64,prevStepVal::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
  exacteA(q,cacheA,index,index)
  aii=cacheA[1]
  exacteA(q,cacheA,j,j)
  ajj=cacheA[1]
  exacteA(q,cacheA,index,j)
  aij=cacheA[1]
  exacteA(q,cacheA,j,index)
  aji=cacheA[1]

  uii=dxaux[index][1]-aii*qaux[index][1]
 
  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;
 
 
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;
  xj=x[j][0]
  tx[j]=simt
 
      
     
     
     #=  uij=dxaux[index][1]-aii*qaux[index][1]-aij*qj#qaux[j][1]
  uji=ẋj-ajj*qj-aji*qaux[index][1]
  dxj=aji*qi+ajj*qj+uji  =#
  ujj=ẋj-ajj*qj
  #u[j][j][1]=ujj
  uji=ujj-aji*qaux[index][1]# 
  dxj=aji*qi+ajj*qaux[j][1]+uji
  iscycle=false 
  if  ((abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj))#= && abs(ẋj * dxj) > 1e-3  =# )|| (dxj*ẋj)<=0.0
  #if (aji*(qi-qaux[index][1])+ẋj)*ẋj<0.0
  uij=uii-aij*qaux[j][1]
    qjplus=xj+sign(ẋj)*quanj
    dxi=aii*qi+aij*qjplus+uij
  # if dxi*ẋi<0.0
  # if ẋi*(ẋi+aii*(qi-qaux[index][1])+aij*(xj-qj-sign(ẋj)*quanj))<0.0
  if ((abs(dxi)*3<abs(ẋi) || abs(dxi)>3*abs(ẋi))#= && abs(ẋi * dxi)> 1e-3 =#  )|| (dxi*ẋi)<=0.0
 # if ẋi*(ẋi+aii*(qi-qaux[index][1])+aij*(-2*sign(ẋj)*quanj))<0.0
      iscycle=true  
      
      
   


      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
      qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) 
        h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > 2*quani || abs(qj - xj) > 2*quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (1.96*quani / abs(qi - xi));
        h2 = h * (1.96*quanj / abs(qj - xj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          #println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        if maxIter < 1
          println("maxiter of updateQ      = ",maxIter)
         end
        end

      q[index][0]=qi# store back helper vars
      q[j][0]=qj
      tq[j]=simt

    end #end second dependecy check
 end # end outer dependency check
 return iscycle
end  =#

#= function nmisCycle_and_simulUpdate(::Val{1},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
  exacteA(q,cacheA,index,index)
  aii=cacheA[1]
  exacteA(q,cacheA,j,j)
  ajj=cacheA[1]
  exacteA(q,cacheA,index,j)
  aij=cacheA[1]
  exacteA(q,cacheA,j,index)
  aji=cacheA[1]

  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;#olddx[j][1]=ẋj
     
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;
  xj=x[j][0]
  tx[j]=simt

 
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qaux[index][1]
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qaux[index][1]-aij*qj#qaux[j][1]
    iscycle=false
    dxj=aji*qi+ajj*qaux[j][1]+uji #only future qi
    qjplus=xj+sign(dxj)*quanj

      dxi=aii*qi+aij*qjplus+uij #both future qi & qj
      dxi2=aii*qi+aij*qj+uij #only future qi
      dxj2=aji*qi+ajj*qjplus+uji#both future qi & qj
     # if (abs(ẋj+dxj)<1e-2 ||  dxj2<1e-2 || abs(ẋj)<1e-2 ) && (abs(ẋi+dxi)<1e-2 || abs(dxi2)<1e-2 || abs(ẋi)<1e-2 )
    # if (abs(ẋj*dxj)<0 || (abs(dxj+ẋj)<1e-1  )) && ((ẋi*dxi)<0 || (abs(dxi+ẋi)<1e-1 ))
        
    #  if (abs(ẋj*dxj)<0 ) && ((ẋi*dxi)<0)
        #= dxj=aji*qi+ajj*qaux[j][1]+uji
    iscycle=false
        
    if (dxj*ẋj)<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
     # uij=uii-aij*qaux[j][1]
      qjplus=xj+sign(dxj)*quanj
      dxi=aii*qi+aij*qjplus+uij
      if (dxi*ẋi)<0.0
    =#
   
    #if abs(a[j][index]*(2*quani))>abs(x[j][1])
    #if abs(a[j][index]*(q[index][0]-qaux[index][1]))>abs(x[j][1])
    #if (a[j][index]*(q[index][0]-qaux[index][1])+x[j][1])*x[j][1]<0.0
    #if (dxj*ẋj)<=0.0 || abs(dxj-ẋj)>abs(dxj+ẋj)/2
   #  if  abs(dxj-ẋj)>abs(dxj+ẋj)/2
     if  #= abs(dxj-ẋj)>abs(dxj+ẋj)/1.8 =#((abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj))#= && abs(ẋj * dxj) > 1e-3  =# )|| (dxj*ẋj)<=0.0#= && abs(dxj-ẋj)>abs(dxj+ẋj)/20 =#
    #  if (dxj*ẋj)<=0.0
      # if abs(a[j][index]*2*quantum[index])>abs(x[j][1])
    ############################################################################
   
      #if x[index][1]*(x[index][1]+a[index][index]*(q[index][0]-qaux[index][1])+a[index][j]*(xj-qj-sign(x[j][1])*quantum[j]))<0.0
        #if x[index][1]*(x[index][1]+a[index][index]*(q[index][0]-qaux[index][1])+a[index][j]*(-2*sign(x[j][1])*quantum[j]))<0.0
     # if x[index][1]*(x[index][1]+a[index][j]*(xj-qj-sign(x[j][1])*quantum[j]))<0.0
     #   if x[index][1]*(x[index][1]+a[index][j]*(-2*sign(x[j][1])*quantum[j]))<0.0
      #= qjplus=xj+sign(dxj)*quanj
      dxi=aii*qi+aij*qjplus+uij =#
     # if (dxi*ẋi)<=0.0 || abs(dxi-ẋi)>abs(dxi+ẋi)/2
    # if  abs(dxi-ẋi)>abs(dxi+ẋi)/2
     if #= abs(dxi-ẋi)>abs(dxi+ẋi)/1.8  =#((abs(dxi)*3<abs(ẋi) || abs(dxi)>3*abs(ẋi))#= && abs(ẋi * dxi)> 1e-3 =#  )|| (dxi*ẋi)<=0.0 #= && abs(dxi-ẋi)>abs(dxi+ẋi)/20 =#
    #if (dxi*ẋi)<=0.0 
    iscycle=true
       #=  qi=xi
        qj=xj =#
        h_two=-1.0;h_three=-1.0
        #= qj=(-aii*uji+aji*uij)/(ajj*aii-aji*aij)
        qi=(-aij*(aji*uij-uji*aii))/(aii*(ajj*aii-aji*aij))-uij/aii =#
       #=  h = (100.0+ft-simt)
        h_one=h
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
        qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        if (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj)  =#
          #h = ft-simt
          h=firstguessH
          h_one=h
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
          qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        #end
      #=   @show h,index,qi,xi
        @show qj,xj
        htest=100.0-simt
        qitest = ((1-htest*ajj)*(xi+htest*uij)+htest*aij*(xj+htest*uji))/Δ
        qjtest = ((1-htest*aii)*(xj+htest*uji)+htest*aji*(xi+htest*uij))/Δ
        @show qitest,qjtest =#
        if (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) #checking qi-xi is not needed since firstguess just made it less than delta
          h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
          h=min(h1,h2)
          h_two=h
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          if Δ==0
            Δ=1e-12
          end
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
          qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        end
        maxIter=1000
        while (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) && (maxIter>0)
            maxIter-=1
            h1 = h * (0.99*quani / abs(qi - xi));
            h2 = h * (0.99*quanj / abs(qj - xj));
            h=min(h1,h2)
            h_three=h
            Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
            if Δ==0
              Δ=1e-12
              println("delta liqss1 simulupdate==0")
            end
            qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
            qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
            if maxIter < 1 println("maxiter of updateQ      = ",maxIter) end
        end 
      #=   if h==h_one
          dxi=aii*qi+aij*qj+uij
          dxj=aji*qi+ajj*qj+uji
          @show index,h,simt,dxi,dxj
        
        elseif h==h_two
          dxi=aii*qi+aij*qj+uij
          dxj=aji*qi+ajj*qj+uji
          @show index,simt,h,quani,quanj,ẋi,ẋj,dxi,dxj
        
          
        else
          dxi=aii*qi+aij*qj+uij
          dxj=aji*qi+ajj*qj+uji
          @show h,index,simt ,abs(qi - xi) , abs(qj - xj),dxi,dxj
          
        end =#
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        tq[j]=simt 
      end #end second dependecy check
   end # end outer dependency check
   return iscycle
end =#

#= function nmisCycle_and_simulUpdate(simuldeltaiVals,simuldeltajVals,simulqxiVals,simulqxjVals, simulHTimes,simulHVals,::Val{1},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
  exacteA(q,cacheA,index,index)
  aii=cacheA[1]
  exacteA(q,cacheA,j,j)
  ajj=cacheA[1]
  exacteA(q,cacheA,index,j)
  aij=cacheA[1]
  exacteA(q,cacheA,j,index)
  aji=cacheA[1]

  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;#olddx[j][1]=ẋj
     
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;
  xj=x[j][0]
  tx[j]=simt

 
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qaux[index][1]
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qaux[index][1]-aij*qj#qaux[j][1]
    iscycle=false
    dxj=aji*qi+ajj*qaux[j][1]+uji #only future qi
    qjplus=xj+sign(dxj)*quanj

      dxi=aii*qi+aij*qjplus+uij #both future qi & qj
      dxi2=aii*qi+aij*qj+uij #only future qi
      dxj2=aji*qi+ajj*qjplus+uji#both future qi & qj
 
        
   
     if  #= abs(dxj-ẋj)>abs(dxj+ẋj)/1.8 =#((abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj))#= && abs(ẋj * dxj) > 1e-3  =# )|| (dxj*ẋj)<=0.0#= && abs(dxj-ẋj)>abs(dxj+ẋj)/20 =#
    #  if (dxj*ẋj)<=0.0
      # if abs(a[j][index]*2*quantum[index])>abs(x[j][1])
    ############################################################################
   
     if #= abs(dxi-ẋi)>abs(dxi+ẋi)/1.8  =#((abs(dxi)*3<abs(ẋi) || abs(dxi)>3*abs(ẋi))#= && abs(ẋi * dxi)> 1e-3 =#  )|| (dxi*ẋi)<=0.0 #= && abs(dxi-ẋi)>abs(dxi+ẋi)/20 =#
    #if (dxi*ẋi)<=0.0 
    iscycle=true
      
        h_two=-1.0;h_three=-1.0
        h=firstguessH
         #=  
          h_one=h
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
          qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ =#

           #    if (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj)
                                      ppi=pointer(Vector{NTuple{2,Float64}}(undef, 7))
                                      ppj=pointer(Vector{NTuple{2,Float64}}(undef, 7))
                                      resppi = pointer(Vector{Float64}(undef, 2))
                                      resppj = pointer(Vector{Float64}(undef, 2))
                                      resppi2 = pointer(Vector{Float64}(undef, 2))
                                      resppj2 = pointer(Vector{Float64}(undef, 2))
                                      unsafe_store!(resppi, -1.0, 1);unsafe_store!(resppi, -1.0, 2)
                                      unsafe_store!(resppi2, -1.0, 1);unsafe_store!(resppi2, -1.0, 2)
                                      unsafe_store!(resppj, -1.0, 1);unsafe_store!(resppj, -1.0, 2)
                                      unsafe_store!(resppj2, -1.0, 1);unsafe_store!(resppj2, -1.0, 2)
                                      #find positive zeros f=+-Δ
                                      bi=aii*xi+aij*xj+uij;ci=aij*(aji*xi+uji)-ajj*(aii*xi+uij);αi=-ajj-aii;βi=aii*ajj-aij*aji
                                      #= coefi=@SVector [quani,αi*quani-bi,βi*quani-ci]#
                                      posSolPlusi= PosRoot(coefi, Val(2)) =#
                                      coefi=NTuple{3,Float64}((βi*quani-ci,αi*quani-bi,quani))
                                      allrealrootintervalnewtonregulafalsi(coefi,resppi,ppi)
                                      resTupi=(unsafe_load(resppi,1),unsafe_load(resppi,2))
                                      posSolPlusi=filter((x) -> x >0.0 , resTupi)
                                    #=  coefi=@SVector [-quani,-αi*quani-bi,-βi*quani-ci]
                                      posSolminusi= PosRoot(coefi, Val(2)) =#
                                      coefi2=NTuple{3,Float64}((-βi*quani-ci,-αi*quani-bi,-quani))
                                      allrealrootintervalnewtonregulafalsi(coefi2,resppi2,ppi)
                                      resTupi2=(unsafe_load(resppi2,1),unsafe_load(resppi2,2))
                                      posSolminusi=filter((x) -> x >0.0 , resTupi2)
                                      
                                      posSoli=(posSolPlusi...,posSolminusi...)

                                      bj=ajj*xj+aji*xi+uji;cj=aji*(aij*xj+uij)-aii*(ajj*xj+uji);αj=-aii-ajj;βj=ajj*aii-aji*aij

                                      #= coefj=@SVector [quanj,αj*quanj-bj,βj*quanj-cj]#
                                      posSolPlusj= PosRoot(coefj, Val(2)) =#
                                      coefj=NTuple{3,Float64}((βj*quanj-cj,αj*quanj-bj,quanj))
                                      allrealrootintervalnewtonregulafalsi(coefj,resppj,ppj)
                                      resTupj=(unsafe_load(resppj,1),unsafe_load(resppj,2))
                                      posSolPlusj=filter((x) -> x >0.0 , resTupj)
                                      #= coefj2=@SVector [-quanj,-αj*quanj-bj,-βj*quanj-cj]#
                                      posSolminusj= PosRoot(coefj2, Val(2)) =#
                                      coefj2=NTuple{3,Float64}((-βj*quanj-cj,-αj*quanj-bj,-quanj))
                                      allrealrootintervalnewtonregulafalsi(coefj2,resppj2,ppj)
                                      resTupj2=(unsafe_load(resppj2,1),unsafe_load(resppj2,2))
                                      posSolminusj=filter((x) -> x >0.0 , resTupj2)

                                      posSolj=(posSolPlusj...,posSolminusj...)

                                      #construct intervals
                                      resi=constructIntrval(posSoli)
                                      resj=constructIntrval(posSolj)
                                   #=  if simt < 0.0 
                                      @show coefj,coefj2
                                      @show posSolPlusj
                                      @show posSolminusj
                                      @show posSoli,posSolj
                                      @show resi,resj
                                    end =#
                                    
                                      #find best H (largest overlap)
                                      while true
                                        if resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
                                        
                                          h=min(resj[end][2],resi[end][2] )
                                          if h==Inf && (resj[end][1]!=0.0 || resi[end][1]!=0.0) # except case both  0 sols  h1=(0,Inf)
                                            h=max(resj[end][1],resi[end][1],ft-simt,firstguessH) #ft-simt in case they are both ver small...elaborate on his later
                                          end
                                          if h==Inf # both lower bounds ==0  --> zero sols for both
                                            #h=ft-simt
                                            qi=xi+ci/βi
                                            qj=xj+cj/βj
                                          else
                                              Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                                              qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                                              qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
                                          end


                                     #=      if abs(qi - xi) / quani< 0.5 && abs(qj - xj) / quanj <0.5 && -0.02999999<simt<-0.028
                                            fhi(y)=(y*bi+y*y*ci)/(1+αi*y+βi*y*y)
                                            fdeltai(y)=-quani
                                            p1i=plot()
                                            p1i=plot!(p1i,fhi)
                                            p1i=plot!(p1i,fdeltai)
                                            savefig(p1i, "plot_simul_adr_FHI_analV0301.png")

                                            
                                            fhj(y)=(y*bj+y*y*cj)/(1+αj*y+βj*y*y)
                                            fdeltaj(y)=-quanj
                                            p1j=plot()
                                            p1j=plot!(fhj)
                                            p1j=plot!(p1j,fdeltaj)
                                            savefig(p1j, "plot_simul_adr_FHJ_analV0301.png")

                                            
                                            #p2j=plot()
                                            p2i=plot(fhi,xlims=(-0.00001,0.5),ylims=(-0.002,-0.0018))
                                            savefig(p2i, "plot_simul_adr_FHi_analV0301_zoomed.png")

                                     
                                            p2j=plot(fhj,xlims=(-0.00001,0.5),ylims=(-0.0073,-0.0068))
                                            savefig(p2j, "plot_simul_adr_FHJ_analV0301_zoomed.png")


                                            @show simt,h
                                            @show qi, xi,quani,qj , xj,quanj
                                            @show abs(qi - xi) / quani, abs(qj - xj) / quanj
                                            @show resi,resj
                                            @show posSoli,posSolj
                                            @show  coefi, coefi2, coefj, coefj2
                                            @show resTupi,resTupi2
                                          end =#
                                        # @show resi,resj
                                          break
                                        else
                                          if resj[end][1]<resi[end][1]#remove last seg of resi
                                            resi=resi[1:end-1]
                                          else#remove last seg of resj
                                            resj=resj[1:end-1]
                                          end
                                        end

                                      end
                                    #=   if h==Inf #h1 is the culprit
                                        #= if ci/βi>0
                                          mi=(quani+ci/βi)/2
                                          coefi=@SVector [mi,αi*mi-bi,βi*mi-ci]# =#
                                          
                                          mi=(ci/βi)/2
                                          coefi=@SVector [mi,αi*mi-bi,βi*mi-ci]
                                          h1i=minPosRoot(coefi, Val(2))
                                          mj=(cj/βj)/2
                                          coefj=@SVector [mj,αj*mj-bj,βj*mj-cj]
                                          h1j=minPosRoot(coefj, Val(2))
                                          h=min(h1i,h1j)
                                         # @show h,h1i,h1j 

                                      end =#
                                     #=  Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                                      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                                      qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ =#



                                      if ((abs(qi - xi) > 1.1*quani || abs(qj - xj) > 1.1*quanj)) #&& simt==4.848661166194424

                                        @show simt,h
                                        @show qi, xi,quani,qj , xj,quanj
                                        @show abs(qi - xi) / quani, abs(qj - xj) / quanj
                                        @show resi,resj
                                        @show posSoli,posSolj
                                        @show  coefi, coefi2, coefj, coefj2
                                    end
             #   end
              #= if simt < 0.0
                @show simt
                @show abs(qi - xi) < quani, abs(qj - xj) < quanj
                @show abs(qi - xi) / quani, abs(qj - xj) / quanj
                @show qi, xi,quani,qj , xj,quanj
              end =#

       
    
       #=  if (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) #checking qi-xi is not needed since firstguess just made it less than delta
          h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
          h=min(h1,h2)
          h_two=h
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          if Δ==0
            Δ=1e-12
          end
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
          qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        end
        maxIter=1000
        while (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) && (maxIter>0)
            maxIter-=1
            h1 = h * (0.99*quani / abs(qi - xi));
            h2 = h * (0.99*quanj / abs(qj - xj));
            h=min(h1,h2)
            h_three=h
            Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
            if Δ==0
              Δ=1e-12
              println("delta liqss1 simulupdate==0")
            end
            qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
            qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
            if maxIter < 1 println("maxiter of updateQ      = ",maxIter) end
        end  =#
       # if maxIter < 950 println("maxiter      = ",maxIter) end
       # @show simt,index,h
         push!(simulHTimes,simt)
       push!(simulHVals,h)
   
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        push!(simulqxiVals,abs(qi-xi))
        push!(simulqxjVals,abs(qj-xj))
        push!(simuldeltaiVals,quani)
        push!(simuldeltajVals,quanj)
        tq[j]=simt 
      end #end second dependecy check
   end # end outer dependency check
   return iscycle
end =#


function nmisCycle_and_simulUpdate(simuldeltaiVals,simuldeltajVals,simulqxiVals,simulqxjVals, simulHTimes,simulHVals,::Val{1},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
  exacteA(q,cacheA,index,index)
  aii=cacheA[1]
  exacteA(q,cacheA,j,j)
  ajj=cacheA[1]
  exacteA(q,cacheA,index,j)
  aij=cacheA[1]
  exacteA(q,cacheA,j,index)
  aji=cacheA[1]

  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;#olddx[j][1]=ẋj
     
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;
  xj=x[j][0]
  tx[j]=simt

 
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qaux[index][1]
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qaux[index][1]-aij*qj#qaux[j][1]
    iscycle=false
    dxj=aji*qi+ajj*qaux[j][1]+uji #only future qi
    qjplus=xj+sign(dxj)*quanj

      dxi=aii*qi+aij*qjplus+uij #both future qi & qj
      dxi2=aii*qi+aij*qj+uij #only future qi
      dxj2=aji*qi+ajj*qjplus+uji#both future qi & qj
 
        
   
     if  #= abs(dxj-ẋj)>abs(dxj+ẋj)/1.8 =#((abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj))#= && abs(ẋj * dxj) > 1e-3  =# )|| (dxj*ẋj)<=0.0#= && abs(dxj-ẋj)>abs(dxj+ẋj)/20 =#
    #  if (dxj*ẋj)<=0.0
      # if abs(a[j][index]*2*quantum[index])>abs(x[j][1])
    ############################################################################
   
     if #= abs(dxi-ẋi)>abs(dxi+ẋi)/1.8  =#((abs(dxi)*3<abs(ẋi) || abs(dxi)>3*abs(ẋi))#= && abs(ẋi * dxi)> 1e-3 =#  )|| (dxi*ẋi)<=0.0 #= && abs(dxi-ẋi)>abs(dxi+ẋi)/20 =#
    #if (dxi*ẋi)<=0.0 
    iscycle=true
      
        h_two=-1.0;h_three=-1.0
       #=  h=firstguessH
         #=  
          h_one=h
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
          qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ =#

           #    if (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj)
                                      ppi=pointer(Vector{NTuple{2,Float64}}(undef, 7))
                                      ppj=pointer(Vector{NTuple{2,Float64}}(undef, 7))
                                      resppi = pointer(Vector{Float64}(undef, 2))
                                      resppj = pointer(Vector{Float64}(undef, 2))
                                      resppi2 = pointer(Vector{Float64}(undef, 2))
                                      resppj2 = pointer(Vector{Float64}(undef, 2))
                                      unsafe_store!(resppi, -1.0, 1);unsafe_store!(resppi, -1.0, 2)
                                      unsafe_store!(resppi2, -1.0, 1);unsafe_store!(resppi2, -1.0, 2)
                                      unsafe_store!(resppj, -1.0, 1);unsafe_store!(resppj, -1.0, 2)
                                      unsafe_store!(resppj2, -1.0, 1);unsafe_store!(resppj2, -1.0, 2)
                                      #find positive zeros f=+-Δ
                                      bi=aii*xi+aij*xj+uij;ci=aij*(aji*xi+uji)-ajj*(aii*xi+uij);αi=-ajj-aii;βi=aii*ajj-aij*aji
                                      #= coefi=@SVector [quani,αi*quani-bi,βi*quani-ci]#
                                      posSolPlusi= PosRoot(coefi, Val(2)) =#
                                      coefi=NTuple{3,Float64}((βi*quani-ci,αi*quani-bi,quani))
                                      allrealrootintervalnewtonregulafalsi(coefi,resppi,ppi)
                                      resTupi=(unsafe_load(resppi,1),unsafe_load(resppi,2))
                                      posSolPlusi=filter((x) -> x >0.0 , resTupi)
                                    #=  coefi=@SVector [-quani,-αi*quani-bi,-βi*quani-ci]
                                      posSolminusi= PosRoot(coefi, Val(2)) =#
                                      coefi2=NTuple{3,Float64}((-βi*quani-ci,-αi*quani-bi,-quani))
                                      allrealrootintervalnewtonregulafalsi(coefi2,resppi2,ppi)
                                      resTupi2=(unsafe_load(resppi2,1),unsafe_load(resppi2,2))
                                      posSolminusi=filter((x) -> x >0.0 , resTupi2)
                                      
                                      posSoli=(posSolPlusi...,posSolminusi...)

                                      bj=ajj*xj+aji*xi+uji;cj=aji*(aij*xj+uij)-aii*(ajj*xj+uji);αj=-aii-ajj;βj=ajj*aii-aji*aij

                                      #= coefj=@SVector [quanj,αj*quanj-bj,βj*quanj-cj]#
                                      posSolPlusj= PosRoot(coefj, Val(2)) =#
                                      coefj=NTuple{3,Float64}((βj*quanj-cj,αj*quanj-bj,quanj))
                                      allrealrootintervalnewtonregulafalsi(coefj,resppj,ppj)
                                      resTupj=(unsafe_load(resppj,1),unsafe_load(resppj,2))
                                      posSolPlusj=filter((x) -> x >0.0 , resTupj)
                                      #= coefj2=@SVector [-quanj,-αj*quanj-bj,-βj*quanj-cj]#
                                      posSolminusj= PosRoot(coefj2, Val(2)) =#
                                      coefj2=NTuple{3,Float64}((-βj*quanj-cj,-αj*quanj-bj,-quanj))
                                      allrealrootintervalnewtonregulafalsi(coefj2,resppj2,ppj)
                                      resTupj2=(unsafe_load(resppj2,1),unsafe_load(resppj2,2))
                                      posSolminusj=filter((x) -> x >0.0 , resTupj2)

                                      posSolj=(posSolPlusj...,posSolminusj...)

                                      #construct intervals
                                      resi=constructIntrval(posSoli)
                                      resj=constructIntrval(posSolj)
                                   #=  if simt < 0.0 
                                      @show coefj,coefj2
                                      @show posSolPlusj
                                      @show posSolminusj
                                      @show posSoli,posSolj
                                      @show resi,resj
                                    end =#
                                    
                                      #find best H (largest overlap)
                                      while true
                                        if resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
                                        
                                          h=min(resj[end][2],resi[end][2] )
                                          if h==Inf && (resj[end][1]!=0.0 || resi[end][1]!=0.0) # except case both  0 sols  h1=(0,Inf)
                                            h=max(resj[end][1],resi[end][1],ft-simt,firstguessH) #ft-simt in case they are both ver small...elaborate on his later
                                          end
                                          if h==Inf # both lower bounds ==0  --> zero sols for both
                                            #h=ft-simt
                                            qi=xi+ci/βi
                                            qj=xj+cj/βj
                                          else
                                              Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                                              qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                                              qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
                                          end


                                     #=      if abs(qi - xi) / quani< 0.5 && abs(qj - xj) / quanj <0.5 && -0.02999999<simt<-0.028
                                            fhi(y)=(y*bi+y*y*ci)/(1+αi*y+βi*y*y)
                                            fdeltai(y)=-quani
                                            p1i=plot()
                                            p1i=plot!(p1i,fhi)
                                            p1i=plot!(p1i,fdeltai)
                                            savefig(p1i, "plot_simul_adr_FHI_analV0301.png")

                                            
                                            fhj(y)=(y*bj+y*y*cj)/(1+αj*y+βj*y*y)
                                            fdeltaj(y)=-quanj
                                            p1j=plot()
                                            p1j=plot!(fhj)
                                            p1j=plot!(p1j,fdeltaj)
                                            savefig(p1j, "plot_simul_adr_FHJ_analV0301.png")

                                            
                                            #p2j=plot()
                                            p2i=plot(fhi,xlims=(-0.00001,0.5),ylims=(-0.002,-0.0018))
                                            savefig(p2i, "plot_simul_adr_FHi_analV0301_zoomed.png")

                                     
                                            p2j=plot(fhj,xlims=(-0.00001,0.5),ylims=(-0.0073,-0.0068))
                                            savefig(p2j, "plot_simul_adr_FHJ_analV0301_zoomed.png")


                                            @show simt,h
                                            @show qi, xi,quani,qj , xj,quanj
                                            @show abs(qi - xi) / quani, abs(qj - xj) / quanj
                                            @show resi,resj
                                            @show posSoli,posSolj
                                            @show  coefi, coefi2, coefj, coefj2
                                            @show resTupi,resTupi2
                                          end =#
                                        # @show resi,resj
                                          break
                                        else
                                          if resj[end][1]<resi[end][1]#remove last seg of resi
                                            resi=resi[1:end-1]
                                          else#remove last seg of resj
                                            resj=resj[1:end-1]
                                          end
                                        end

                                      end
                                    #=   if h==Inf #h1 is the culprit
                                        #= if ci/βi>0
                                          mi=(quani+ci/βi)/2
                                          coefi=@SVector [mi,αi*mi-bi,βi*mi-ci]# =#
                                          
                                          mi=(ci/βi)/2
                                          coefi=@SVector [mi,αi*mi-bi,βi*mi-ci]
                                          h1i=minPosRoot(coefi, Val(2))
                                          mj=(cj/βj)/2
                                          coefj=@SVector [mj,αj*mj-bj,βj*mj-cj]
                                          h1j=minPosRoot(coefj, Val(2))
                                          h=min(h1i,h1j)
                                         # @show h,h1i,h1j 

                                      end =#
                                     #=  Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                                      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                                      qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ =#



                                      if ((abs(qi - xi) > 1.1*quani || abs(qj - xj) > 1.1*quanj)) #&& simt==4.848661166194424

                                        @show simt,h
                                        @show qi, xi,quani,qj , xj,quanj
                                        @show abs(qi - xi) / quani, abs(qj - xj) / quanj
                                        @show resi,resj
                                        @show posSoli,posSolj
                                        @show  coefi, coefi2, coefj, coefj2
                                    end
             #   end =#
              #= if simt < 0.0
                @show simt
                @show abs(qi - xi) < quani, abs(qj - xj) < quanj
                @show abs(qi - xi) / quani, abs(qj - xj) / quanj
                @show qi, xi,quani,qj , xj,quanj
              end =#

       
     
          h=ft-simt
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
          qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        if (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) #checking qi-xi is not needed since firstguess just made it less than delta
          h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
          h=min(h1,h2)
          h_two=h
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          if Δ==0
            Δ=1e-12
          end
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
          qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
        end
        maxIter=1000
        while (abs(qi - xi) > 1*quani || abs(qj - xj) > 1*quanj) && (maxIter>0)
            maxIter-=1
            h1 = h * (0.99*quani / abs(qi - xi));
            h2 = h * (0.99*quanj / abs(qj - xj));
            h=min(h1,h2)
            h_three=h
            Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
            if Δ==0
              Δ=1e-12
              println("delta liqss1 simulupdate==0")
            end
            qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
            qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
            if maxIter < 1 println("maxiter of updateQ      = ",maxIter) end
        end 
       # if maxIter < 950 println("maxiter      = ",maxIter) end
       # @show simt,index,h
         push!(simulHTimes,simt)
       push!(simulHVals,h)
   
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        push!(simulqxiVals,abs(qi-xi))
        push!(simulqxjVals,abs(qj-xj))
        push!(simuldeltaiVals,quani)
        push!(simuldeltajVals,quanj)
        tq[j]=simt 
      end #end second dependecy check
   end # end outer dependency check
   return iscycle
end