
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
 
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;xjaux=x[j][0]
  tx[j]=simt
  ujj=ẋj-ajj*qj
   
  uji=ujj-aji*qaux[index][1]
 
  dxj=aji*qi+ajj*qj+uji
  iscycle=false 
  if dxj*ẋj<0
    qjplus=xjaux+sign(dxj)*quanj
    uij=uii-aij*qj#qaux[j][1]
   
    dxi=aii*qi+aij*qjplus+uij
    if dxi*ẋi<0
      iscycle=true  
        
      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
      qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
        h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (1.96*quani / abs(qi - xi));
        h2 = h * (1.96*quanj / abs(qj - xjaux));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          #println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
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
  xjaux=x[j][0]
  tx[j]=simt
  #dxj=aji*qi+ajj*qj+uji 
  iscycle=false 
  if (aji*(qi-qaux[index][1])+ẋj)*ẋj<0.0
    exacteA(q,cacheA,index,j)
    aij=cacheA[1]
    exacteA(q,cacheA,index,index)
    aii=cacheA[1]
   # qjplus=xjaux-sign(ẋj)*quanj
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
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
      qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
        h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (1.96*quani / abs(qi - xi));
        h2 = h * (1.96*quanj / abs(qj - xjaux));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          #println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
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
  xjaux=x[j][0]
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
   
    qjplus=xjaux+sign(dxj)*quanj
    dxi=aii*qi+aij*qjplus+uij
  if (dxi*ẋi)<0.0
  #if ẋi*(ẋi+aii*(qi-qaux[index][1])+aij*(-2*sign(ẋj)*quanj))<0.0
      iscycle=true  
      
   
   


      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
      qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
        h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (1.96*quani / abs(qi - xi));
        h2 = h * (1.96*quanj / abs(qj - xjaux));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          #println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
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
  xjaux=x[j][0]
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
    qjplus=xjaux+sign(ẋj)*quanj
    dxi=aii*qi+aij*qjplus+uij
  # if dxi*ẋi<0.0
  # if ẋi*(ẋi+aii*(qi-qaux[index][1])+aij*(xj-qj-sign(ẋj)*quanj))<0.0
  if ((abs(dxi)*3<abs(ẋi) || abs(dxi)>3*abs(ẋi))#= && abs(ẋi * dxi)> 1e-3 =#  )|| (dxi*ẋi)<=0.0
 # if ẋi*(ẋi+aii*(qi-qaux[index][1])+aij*(-2*sign(ẋj)*quanj))<0.0
      iscycle=true  
      
      
   


      h = ft-simt
      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
      qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
        h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
      end
      maxIter=1000
      while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
        maxIter-=1
        h1 = h * (1.96*quani / abs(qi - xi));
        h2 = h * (1.96*quanj / abs(qj - xjaux));
        h=min(h1,h2)
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        if Δ==0
          Δ=1e-12
          #println("delta liqss1 simulupdate==0")
        end
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
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

function nmisCycle_and_simulUpdate(::Val{1},index::Int,j::Int,dirI::Float64,prevStepVal::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
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
  xjaux=x[j][0]
  tx[j]=simt

 
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qaux[index][1]
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qaux[index][1]-aij*qj#qaux[j][1]
    iscycle=false
    dxj=aji*qi+ajj*qaux[j][1]+uji
    #= dxj=aji*qi+ajj*qaux[j][1]+uji
    iscycle=false

    if (dxj*ẋj)<0.0 #(dqjplus*qj1)<=0.0 with dir is better since when dir =0 we do not enter
     # uij=uii-aij*qaux[j][1]
      qjplus=xjaux+sign(dxj)*quanj
      dxi=aii*qi+aij*qjplus+uij
      if (dxi*ẋi)<0.0
    =#
    #if abs(a[j][index]*(2*quani))>abs(x[j][1])
    #if abs(a[j][index]*(q[index][0]-qaux[index][1]))>abs(x[j][1])
    #if (a[j][index]*(q[index][0]-qaux[index][1])+x[j][1])*x[j][1]<0.0
      if  #= abs(dxj-ẋj)>abs(dxj+ẋj)/1.8 =#((abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj))#= && abs(ẋj * dxj) > 1e-3  =# )|| (dxj*ẋj)<=0.0#= && abs(dxj-ẋj)>abs(dxj+ẋj)/20 =#
    # if abs(a[j][index]*2*quantum[index])>abs(x[j][1])
    ############################################################################
   
      #if x[index][1]*(x[index][1]+a[index][index]*(q[index][0]-qaux[index][1])+a[index][j]*(xjaux-qj-sign(x[j][1])*quantum[j]))<0.0
        #if x[index][1]*(x[index][1]+a[index][index]*(q[index][0]-qaux[index][1])+a[index][j]*(-2*sign(x[j][1])*quantum[j]))<0.0
     # if x[index][1]*(x[index][1]+a[index][j]*(xjaux-qj-sign(x[j][1])*quantum[j]))<0.0
     #   if x[index][1]*(x[index][1]+a[index][j]*(-2*sign(x[j][1])*quantum[j]))<0.0
      qjplus=xjaux+sign(dxj)*quanj
      dxi=aii*qi+aij*qjplus+uij
      if #= abs(dxi-ẋi)>abs(dxi+ẋi)/1.8  =#((abs(dxi)*3<abs(ẋi) || abs(dxi)>3*abs(ẋi))#= && abs(ẋi * dxi)> 1e-3 =#  )|| (dxi*ẋi)<=0.0 #= && abs(dxi-ẋi)>abs(dxi+ẋi)/20 =#
        iscycle=true
        h = ft-simt
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
        if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) 
          h1 = (abs(quani / ẋi));h2 = (abs(quanj / ẋj));
          h=min(h1,h2)
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          if Δ==0
            Δ=1e-12
          end
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
          qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
        end
        maxIter=1000
        while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * (1.96*quani / abs(qi - xi));
          h2 = h * (1.96*quanj / abs(qj - xjaux));
          h=min(h1,h2)
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          if Δ==0
            Δ=1e-12
            #println("delta liqss1 simulupdate==0")
          end
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
          qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
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
end