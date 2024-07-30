#analy


 

 function isCycle_and_simulUpdate(::Val{:SimulIter},::Val{Detection},acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,index::Int,j::Int, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)where{Detection}
 
  cacheA[1]=0.0;exacteA(q,cacheA,index,index)
  aii=cacheA[1]
  cacheA[1]=0.0;exacteA(q,cacheA,j,j)
  ajj=cacheA[1]
  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;
     
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;
  xjaux=x[j][0]
  tx[j]=simt

 
   
 
  
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qaux[index][1]
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qaux[index][1]-aij*qj#qaux[j][1]
    iscycle=false
    dxj=aji*qi+ajj*qaux[j][1]+uji
    dxP=aii*qi+aij*qj+uij
    if abs(dxP)<1e-15 && dxP!=0.0
      dxP=0.0
    end
    qjplus=xjaux+sign(dxj)*quanj
      dxi=aii*qi+aij*qjplus+uij
     
      iscycle=detect(Val(Detection),ẋi,dxP,dxi,ẋj,dxj)
      #= if   (dxj*ẋj)<=0.0#= && abs(dxj-ẋj)>abs(dxj+ẋj)/20 =#

      
      if  (dxi*ẋi)<=0.0  =#
    if iscycle
        #iscycle=true
        h = ft-simt
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
        qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
        if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 1*quanj) 
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
        while (abs(qi - xi) > 1*quani || abs(qj - xjaux) > 1*quanj) && (maxIter>0)
          maxIter-=1
          h1 = h * (0.98*quani / abs(qi - xi));
          h2 = h * (0.98*quanj / abs(qj - xjaux));
          h=min(h1,h2)
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          if Δ==0
            Δ=1e-12
            #println("delta liqss1 simulupdate==0")
          end
          qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
          qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
          
        end 
          if maxIter < 1
            println("maxiter of simul_updateQ      = ",maxIter)
            return false
           end
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        tq[j]=simt 
     # end #end second dependecy check
   end # end outer dependency check
   return iscycle
end


function isCycle_and_simulUpdate(::Val{:SimulAna},::Val{Detection},acceptedi::Vector{Vector{Float64}},acceptedj::Vector{Vector{Float64}},aij::Float64,aji::Float64,index::Int,j::Int, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)where{Detection}
 
  cacheA[1]=0.0;exacteA(q,cacheA,index,index)
  aii=cacheA[1]
  cacheA[1]=0.0;exacteA(q,cacheA,j,j)
  ajj=cacheA[1]
  xi=x[index][0];xj=x[j][0];ẋi=x[index][1];ẋj=x[j][1]
  qi=q[index][0];qj=q[j][0]
  quanj=quantum[j];quani=quantum[index]
  qaux[j][1]=qj;
     
  elapsed = simt - tx[j];x[j][0]= xj+elapsed*ẋj;
  xjaux=x[j][0]
  tx[j]=simt

 
   
 
  
   #ujj=ẋj-ajj*qj
    uji=ẋj-ajj*qj-aji*qaux[index][1]
    #uii=dxaux[index][1]-aii*qaux[index][1]
    uij=dxaux[index][1]-aii*qaux[index][1]-aij*qj#qaux[j][1]
    iscycle=false
    dxj=aji*qi+ajj*qaux[j][1]+uji
    dxP=aii*qi+aij*qj+uij
    if abs(dxP)<1e-15 && dxP!=0.0
      dxP=0.0
    end
    qjplus=xjaux+sign(dxj)*quanj
      dxi=aii*qi+aij*qjplus+uij
     
      iscycle=detect(Val(Detection),ẋi,dxP,dxi,ẋj,dxj)
      #= if   (dxj*ẋj)<=0.0#= && abs(dxj-ẋj)>abs(dxj+ẋj)/20 =#

      
      if  (dxi*ẋi)<=0.0  =#
        if iscycle
          #clear accIntrvals and cache of roots
          for i = 1:3# 3 ord1 ,7 ord2
            acceptedi[i][1] = 0.0
            acceptedi[i][2] = 0.0
            acceptedj[i][1] = 0.0
            acceptedj[i][2] = 0.0
          end
       
          #find positive zeros f=+-Δ
          bi = aii * xi + aij * xj + uij
          ci = aij * (aji * xi + uji) - ajj * (aii * xi + uij)
          αi = -ajj - aii
          βi = aii * ajj - aij * aji
          bj = ajj * xj + aji * xi + uji
          cj = aji * (aij * xj + uij) - aii * (ajj * xj + uji)
          αj = -aii - ajj
          βj = ajj * aii - aji * aij
          coefi = NTuple{3,Float64}((βi * quani - ci, αi * quani - bi, quani))
          coefi2 = NTuple{3,Float64}((-βi * quani - ci, -αi * quani - bi, -quani))
          coefj = NTuple{3,Float64}((βj * quanj - cj, αj * quanj - bj, quanj))
          coefj2 = NTuple{3,Float64}((-βj * quanj - cj, -αj * quanj - bj, -quanj))


          if (abs(coefi[1]) > 1e8abs(coefi[2]) ) || (abs(coefi[1]) < 1e-8abs(coefi[2]) ) || (abs(coefj[1]) > 1e8abs(coefj[2]) ) || (abs(coefj[1]) < 1e-8abs(coefj[2]))
                  h = ft-simt
                  Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                  qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
                  qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
                  if (abs(qi - xi) > 1*quani || abs(qj - xjaux) > 1*quanj) 
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
                  while (abs(qi - xi) > 1*quani || abs(qj - xjaux) > 1*quanj) && (maxIter>0)
                    maxIter-=1
                    h1 = h * (0.98*quani / abs(qi - xi));
                    h2 = h * (0.98*quanj / abs(qj - xjaux));
                    h=min(h1,h2)
                    Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                    if Δ==0
                      Δ=1e-12
                      #println("delta liqss1 simulupdate==0")
                    end
                    qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
                    qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
                    
                  end 
                    if maxIter < 1
                      println("maxiter of simul_updateQ      = ",maxIter)
                      return false
                    end
                  q[index][0]=qi# store back helper vars
                  q[j][0]=qj
                  tq[j]=simt 
                  return true
          end
           

          resi1, resi2 = quadRootv2(coefi)
          resi3, resi4 = quadRootv2(coefi2)
          resj1, resj2 = quadRootv2(coefj)
          resj3, resj4 = quadRootv2(coefj2)
          #construct intervals
          constructIntrval(acceptedi, resi1, resi2, resi3, resi4)
          constructIntrval(acceptedj, resj1, resj2, resj3, resj4)
          #find best H (largest overlap)
          ki = 3
          kj = 3
          while true
            currentLi = acceptedi[ki][1]
            currentHi = acceptedi[ki][2]
            currentLj = acceptedj[kj][1]
            currentHj = acceptedj[kj][2]
            if currentLj <= currentLi < currentHj || currentLi <= currentLj < currentHi#resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
              h = min(currentHi, currentHj)
              if h == Inf && (currentLj != 0.0 || currentLi != 0.0) # except case both  0 sols  h1=(0,Inf)
                h = max(currentLj, currentLi) #ft-simt in case they are both ver small...elaborate on his later
              end
              if h == Inf # both lower bounds ==0  --> zero sols for both
                qi = getQfromAsymptote(simt, xi, βi, ci, αi, bi)
                qj = getQfromAsymptote(simt, xj, βj, cj, αj, bj)
              else
                Δ = (1 - h * aii) * (1 - h * ajj) - h * h * aij * aji
                if Δ != 0.0
                  qi = ((1 - h * ajj) * (xi + h * uij) + h * aij * (xj + h * uji)) / Δ
                  qj = ((1 - h * aii) * (xj + h * uji) + h * aji * (xi + h * uij)) / Δ
                else
                  qi, qj, h, maxIter = IterationH(h, xi, quani, xj, quanj, aii, ajj, aij, aji, uij, uji)
                  if maxIter == 0
                    return false
                  end
                end
              end
              break
            else
              if currentHj == 0.0 && currentHi == 0.0#empty
                ki -= 1
                kj -= 1
              elseif currentHj == 0.0 && currentHi != 0.0
                kj -= 1
              elseif currentHi == 0.0 && currentHj != 0.0
                ki -= 1
              else #both non zero
                if currentLj < currentLi#remove last seg of resi
                  ki -= 1
                else #remove last seg of resj
                  kj -= 1
                end
              end
            end
          end # end while       
          q[j][0] = qj
          tq[j] = simt
          q[index][0] = qi# store back helper vars
       
        end
        return iscycle
end


  
function getQfromAsymptote(simt, x::Float64, β::P, c::P, α::P, b::P) where {P<:Union{BigFloat,Float64}} 
  q = 0.0
  if β == 0.0 && c != 0.0
    println("report bug: β==0 && c!=0.0: this leads to asym=Inf while this code came from asym<Delta at $simt")
  elseif β == 0.0 && c == 0.0
    if α == 0.0 && b != 0.0
      println("report bug: α==0 && b!=0.0 this leads to asym=Inf while this code came from asym<Delta at $simt")
    elseif α == 0.0 && b == 0.0
      q = x
    else
      q = x + b / α
    end
  else
    q = x + c / β #asymptote...even though not guaranteed to be best
  end
  q
end
function iterationH(h::Float64, xi::Float64, quani::Float64, xj::Float64, quanj::Float64, aii::Float64, ajj::Float64, aij::Float64, aji::Float64, uij::Float64, uji::Float64)
  qi, qj = 0.0, 0.0
  maxIter = 1000
  while (abs(qi - xi) > 1 * quani || abs(qj - xj) > 1 * quanj) && (maxIter > 0)
    maxIter -= 1
    h1 = h * (0.99 * quani / abs(qi - xi))
    h2 = h * (0.99 * quanj / abs(qj - xj))
    h = min(h1, h2)
    h_three = h
    Δ = (1 - h * aii) * (1 - h * ajj) - h * h * aij * aji
    if Δ != 0
      qi = ((1 - h * ajj) * (xi + h * uij) + h * aij * (xj + h * uji)) / Δ
      qj = ((1 - h * aii) * (xj + h * uji) + h * aji * (xi + h * uij)) / Δ
    end
    # if maxIter < 1 println("maxiter of updateQ      = ",maxIter) end
  end
  qi, qj, h, maxIter
end 

 function detect(::Val{0},ẋi,dxP,dxi,ẋj,dxj)
 
    return false
 
 end
 function detect(::Val{1},ẋi,dxP,dxi,ẋj,dxj)
  if (dxj*ẋj)<=0.0 && (dxi*ẋi)<=0.0
    
      return true
    
  else  
    return false
  end
 end
 function detect(::Val{2},ẋi,dxP,dxi,ẋj,dxj)
  if (dxj*ẋj)<=0.0 && (dxi*dxP)<=0.0
    
      return true
    
  else  
    #@show dxj,ẋj,dxi,dxP
    return false
  end
 end

 function detect(::Val{3},ẋi,dxP,dxi,ẋj,dxj)
  if abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-ẋi)>abs(dxi+ẋi)/2
   
      return true
   
  else  
    return false
  end
 end

 function detect(::Val{4},ẋi,dxP,dxi,ẋj,dxj)
  if abs(dxj-ẋj)>abs(dxj+ẋj)/2 && abs(dxi-dxP)>abs(dxi+dxP)/2
   
      return true
   
  else  
    return false
  end
 end


 