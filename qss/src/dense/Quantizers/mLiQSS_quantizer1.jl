#analy
function nmisCycle_and_simulUpdate(cacherealPosi::Vector{Vector{Float64}},cacherealPosj::Vector{Vector{Float64}},aij::Float64,aji::Float64,respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},trackSimul,::Val{1},index::Int,j::Int,dirI::Float64,dti::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
 
  exacteA(q,cacheA,index,index)
  aii=cacheA[1]
  exacteA(q,cacheA,j,j)
  ajj=cacheA[1]
 #=  exacteA(q,cacheA,index,j)
  aij=cacheA[1]
  exacteA(q,cacheA,j,index)
  aji=cacheA[1] =#

  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  #qaux[j][1]=qj;
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;

  xj=x[j][0]
  tx[j]=simt

 qiminus=qaux[index][1]
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qiminus
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qiminus-aij*qj
    iscycle=false
    dxj=aji*qi+ajj*qj+uji #only future qi   #emulate fj
   
    

    dxithrow=aii*qi+aij*qj+uij #only future qi
                                                      
  qjplus=xj+sign(dxj)*quanj  #emulate updateQ(j)...

    dxi=aii*qi+aij*qjplus+uij #both future qi & qj   #emulate fi
    #dxj2=ajj*qjplus+aji*qi+uji
    
  
    if abs(dxithrow)<1e-15 && dxithrow!=0.0
      dxithrow=0.0
    end
                             
   ########condition:Union 
#=  if abs(dxj)*3<abs(ẋj) || abs(dxj)>3*abs(ẋj) || (dxj*ẋj)<0.0 
    if abs(dxi)>3*abs(ẋi) || abs(dxi)*3<abs(ẋi) ||  (dxi*ẋi)<0.0 
        iscycle=true
    end
  end    =#                          
    ########condition:Union i
    if abs(dxj-ẋj)>(abs(dxj+ẋj)/2)  
      if abs(dxi-dxithrow)>(abs(dxi+dxithrow)/2) 
        iscycle=true
      end
    end
 

   if iscycle
      #trackSimul[1]+=1 # do not have to recomputeNext
       #find positive zeros f=+-Δ
      bi=aii*xi+aij*xj+uij;ci=aij*(aji*xi+uji)-ajj*(aii*xi+uij);αi=-ajj-aii;βi=aii*ajj-aij*aji
      bj=ajj*xj+aji*xi+uji;cj=aji*(aij*xj+uij)-aii*(ajj*xj+uji);αj=-aii-ajj;βj=ajj*aii-aji*aij
    
      coefi=NTuple{3,Float64}((βi*quani-ci,αi*quani-bi,quani))
      coefi2=NTuple{3,Float64}((-βi*quani-ci,-αi*quani-bi,-quani))
      coefj=NTuple{3,Float64}((βj*quanj-cj,αj*quanj-bj,quanj))
      coefj2=NTuple{3,Float64}((-βj*quanj-cj,-αj*quanj-bj,-quanj))

      resi1,resi2=quadRootv2(coefi)
      resi3,resi4=quadRootv2(coefi2)
      resj1,resj2=quadRootv2(coefj)
      resj3,resj4=quadRootv2(coefj2)
     #=  unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2)
      allrealrootintervalnewtonregulafalsi(coefi,respp,pp)
      resi1,resi2=unsafe_load(respp,1),unsafe_load(respp,2) 
   
      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2)
      allrealrootintervalnewtonregulafalsi(coefi2,respp,pp)
      resi3,resi4=unsafe_load(respp,1),unsafe_load(respp,2) 
  
      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2)
      allrealrootintervalnewtonregulafalsi(coefj,respp,pp)
      resj1,resj2=unsafe_load(respp,1),unsafe_load(respp,2) 
     
      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2)
      allrealrootintervalnewtonregulafalsi(coefj2,respp,pp)
      resj3,resj4=unsafe_load(respp,1),unsafe_load(respp,2)  =#
      



      #construct intervals
      constructIntrval(cacherealPosi,resi1,resi2,resi3,resi4)

      constructIntrval(cacherealPosj,resj1,resj2,resj3,resj4)

      #find best H (largest overlap)
      ki=3;kj=3
      while true
      currentLi=cacherealPosi[ki][1];currentHi=cacherealPosi[ki][2]
      currentLj=cacherealPosj[kj][1];currentHj=cacherealPosj[kj][2]
      if currentLj<=currentLi<currentHj || currentLi<=currentLj<currentHi#resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
                  h=min(currentHi,currentHj )
                  if h==Inf && (currentLj!=0.0 || currentLi!=0.0) # except case both  0 sols  h1=(0,Inf)
                    h=max(currentLj,currentLi#= ,ft-simt,dti =#) #ft-simt in case they are both ver small...elaborate on his later
                  end
                  if h==Inf # both lower bounds ==0  --> zero sols for both
                    qi=xi+ci/βi
                    qj=xj+cj/βj
                  else
                      Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                      qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xj+h*uji))/Δ
                      qj = ((1-h*aii)*(xj+h*uji)+h*aji*(xi+h*uij))/Δ
                  end
                  break
      else
                if currentHj==0.0 && currentHi==0.0#empty
                  ki-=1;kj-=1
                elseif currentHj==0.0 && currentHi!=0.0
                  kj-=1
                elseif currentHi==0.0 && currentHj!=0.0
                  ki-=1
                else #both non zero
                      if currentLj<currentLi#remove last seg of resi
                        ki-=1
                      else #remove last seg of resj
                        kj-=1
                      end
                end
      end

      end
          
      q[j][0]=qj
      q[index][0]=qi# store back helper vars
      tq[j]=simt 
   end
  return iscycle
end  


