
function nmisCycle_and_simulUpdate(cacherealPosi::Vector{Vector{Float64}},cacherealPosj::Vector{Vector{Float64}},respp::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},gh,#= simuldeltaiVals,simuldeltajVals,simulqxiVals,simulqxjVals, simulHTimes,simulHVals, =#::Val{1},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
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
                                     #=  pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
                                    
                                      respp = pointer(Vector{Float64}(undef, 2)) =#
                                      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2) #force negative values
                                      
                                      #find positive zeros f=+-Δ
                                      bi=aii*xi+aij*xj+uij;ci=aij*(aji*xi+uji)-ajj*(aii*xi+uij);αi=-ajj-aii;βi=aii*ajj-aij*aji
                                      #= coefi=@SVector [quani,αi*quani-bi,βi*quani-ci]#
                                      posSolPlusi= PosRoot(coefi, Val(2)) =#
                                      coefi=NTuple{3,Float64}((βi*quani-ci,αi*quani-bi,quani))
                                      allrealrootintervalnewtonregulafalsi(coefi,respp,pp)
                                      #= resTupi=(unsafe_load(respp,1),unsafe_load(respp,2))
                                      posSolPlusi=filter((x) -> x >0.0 , resTupi) =#
                                      resi1=unsafe_load(respp,1)
                                      resi2=unsafe_load(respp,2)
 
                                      coefi2=NTuple{3,Float64}((-βi*quani-ci,-αi*quani-bi,-quani))
                                      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2) #force negative values
                                      allrealrootintervalnewtonregulafalsi(coefi2,respp,pp)
                                  
                                      resi3=unsafe_load(respp,1)
                                      resi4=unsafe_load(respp,2)
                                   

                                      bj=ajj*xj+aji*xi+uji;cj=aji*(aij*xj+uij)-aii*(ajj*xj+uji);αj=-aii-ajj;βj=ajj*aii-aji*aij

                                    
                                      coefj=NTuple{3,Float64}((βj*quanj-cj,αj*quanj-bj,quanj))
                                       unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2) #force negative values
                                       allrealrootintervalnewtonregulafalsi(coefj,respp,pp)
                                       resj1=unsafe_load(respp,1)
                                       resj2=unsafe_load(respp,2)
                                   
                                       coefj2=NTuple{3,Float64}((-βj*quanj-cj,-αj*quanj-bj,-quanj))
                                      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2) #force negative values
                                      allrealrootintervalnewtonregulafalsi(coefj2,respp,pp)
                                      resj3=unsafe_load(respp,1)
                                      resj4=unsafe_load(respp,2)
                                  
                                      #construct intervals
                                      constructIntrval(cacherealPosi,resi1,resi2,resi3,resi4)

                                      constructIntrval(cacherealPosj,resj1,resj2,resj3,resj4)
                                
                                      #find best H (largest overlap)
                            ki=3;kj=3
                          #  @show cacherealPosi
                           # @show cacherealPosj
                          while true
                             # @show ki,kj
                              currentLi=cacherealPosi[ki][1];currentHi=cacherealPosi[ki][2]
                              currentLj=cacherealPosj[kj][1];currentHj=cacherealPosj[kj][2]
                            #  @show currentLi,currentLj,currentHi,currentHj
                              if currentLj<=currentLi<currentHj || currentLi<=currentLj<currentHi#resj[end][1]<=resi[end][1]<=resj[end][2] || resi[end][1]<=resj[end][1]<=resi[end][2] #overlap
                              
                                            h=min(currentHi,currentHj )
                                         #  @show h, currentLi,currentLj,currentHi,currentHj
                                            if h==Inf && (currentLj!=0.0 || currentLi!=0.0) # except case both  0 sols  h1=(0,Inf)
                                              h=max(currentLj,currentLi,ft-simt,firstguessH) #ft-simt in case they are both ver small...elaborate on his later
                                             # @show simt,h
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
                                      



                                    
            
       
    
                                   
   
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
      
        tq[j]=simt 
      end #end second dependecy check
   end # end outer dependency check
   return iscycle
end  

#=  function nmisCycle_and_simulUpdate(gh,#= simuldeltaiVals,simuldeltajVals,simulqxiVals,simulqxjVals, simulHTimes,simulHVals, =#::Val{1},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
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
                                      unsafe_store!(resppi, -1.0, 1);unsafe_store!(resppi, -1.0, 2) #force negative values
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


                                          break
                                        else
                                          if resj[end][1]<resi[end][1]#remove last seg of resi
                                            resi=resi[1:end-1]
                                          else#remove last seg of resj
                                            resj=resj[1:end-1]
                                          end
                                        end

                                      end
                                      



                                    
            
       
    
     
   
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
      
        tq[j]=simt 
      end #end second dependecy check
   end # end outer dependency check
   return iscycle
end   =# 


#=   function nmisCycle_and_simulUpdate(gh,#= simuldeltaiVals,simuldeltajVals,simulqxiVals,simulqxjVals, simulHTimes,simulHVals, =#::Val{1},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
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
      
     
          h=ft-simt
         # h=firstguessH
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
           #=  Δtemp=(1-h1*aii)*(1-h1*ajj)-h1*h1*aij*aji
            qitemp = ((1-h1*ajj)*(xi+h1*uij)+h1*aij*(xj+h1*uji))/Δtemp =#
            h2 = h * (0.99*quanj / abs(qj - xj));
            #= Δtemp=(1-h2*aii)*(1-h2*ajj)-h2*h2*aij*aji
            qjtemp = ((1-h2*ajj)*(xi+h2*uij)+h2*aij*(xj+h2*uji))/Δtemp
            if abs(qitemp - xi) > 1*quani || abs(qjtemp - xj) > 1*quanj
              println("regle de croix did not work")
            end =#
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
         #= push!(simulHTimes,simt)
       push!(simulHVals,h) =#
   
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        #= push!(simulqxiVals,abs(qi-xi))
        push!(simulqxjVals,abs(qj-xj))
        push!(simuldeltaiVals,quani)
        push!(simuldeltajVals,quanj) =#
        tq[j]=simt 
      end #end second dependecy check
   end # end outer dependency check
   return iscycle
end   =#



#= 
function nmisCycle_and_simulUpdate(temporaryhelper,#= simuldeltaiVals,simuldeltajVals,simulqxiVals,simulqxjVals, simulHTimes,simulHVals, =#::Val{1},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64) 
    
    
    
    
  
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
           bi=aii*xi+aij*xj+uij;ci=aij*(aji*xi+uji)-ajj*(aii*xi+uij);αi=-ajj-aii;βi=aii*ajj-aij*aji
           bj=ajj*xj+aji*xi+uji;cj=aji*(aij*xj+uij)-aii*(ajj*xj+uji);αj=-aii-ajj;βj=ajj*aii-aji*aij
           ppi=pointer(Vector{NTuple{2,Float64}}(undef, 7))
           ppj=pointer(Vector{NTuple{2,Float64}}(undef, 7))
           resppi = pointer(Vector{Float64}(undef, 2))
           resppj = pointer(Vector{Float64}(undef, 2))
           resppi2 = pointer(Vector{Float64}(undef, 2))
           resppj2 = pointer(Vector{Float64}(undef, 2))
           unsafe_store!(resppi, -1.0, 1);unsafe_store!(resppi, -1.0, 2) #force negative values
           unsafe_store!(resppi2, -1.0, 1);unsafe_store!(resppi2, -1.0, 2)
           unsafe_store!(resppj, -1.0, 1);unsafe_store!(resppj, -1.0, 2)
           unsafe_store!(resppj2, -1.0, 1);unsafe_store!(resppj2, -1.0, 2)
                 if abs(ci/βi)<quani && abs(cj/βj)<quanj && abs(ci/βi)>quani*0.8 && abs(cj/βj)>quanj*0.8
                  #println("easystep")
                  temporaryhelper[1]+=1
                 
                  qi=xi+ci/βi;qj=xj+cj/βj
                                                       #=  if abs(ci/βi)<quani/2 
                                                                if ci/βi>0
                                                                  coefi=NTuple{3,Float64}((βi*quani-ci,αi*quani-bi,quani))
                                                                  allrealrootintervalnewtonregulafalsi(coefi,resppi,ppi)
                                                                  resTupi=(unsafe_load(resppi,1),unsafe_load(resppi,2))
                                                                  posSolPlusi=filter((x) -> x >0.0 , resTupi)
                                                                          if length(posSolPlusi)>0
                                                                            hi=(max(posSolPlus...,ft-simt,firstguessH),Inf)
                                                                          else
                                                                            hi=(0,Inf)
                                                                          end
                                                                else
                                                                  coefi=NTuple{3,Float64}((-βi*quani-ci,-αi*quani-bi,-quani))
                                                                  allrealrootintervalnewtonregulafalsi(coefi,resppi,ppi)
                                                                  resTupi=(unsafe_load(resppi,1),unsafe_load(resppi,2))
                                                                  posSolPlusi=filter((x) -> x >0.0 , resTupi)
                                                                        if length(posSolPlusi)>0
                                                                          hi=(max(posSolPlusi...,ft-simt,firstguessH),Inf)
                                                                        else
                                                                          hi=(0,Inf)
                                                                        end
                                                                end
                                                        end
                                                        if abs(cj/βj)<quanj/2 
                                                                    if cj/βj>0
                                                                      coefj=NTuple{3,Float64}((βj*quanj-cj,αj*quanj-bj,quanj))
                                                                      allrealrootintervalnewtonregulafalsi(coefj,resppj,ppj)
                                                                      resTupj=(unsafe_load(resppj,1),unsafe_load(resppj,2))
                                                                      posSolPlusj=filter((x) -> x >0.0 , resTupj)
                                                                            if length(posSolPlusj)>0
                                                                              hj=(max(posSolPlusj...,ft-simt,firstguessH),Inf)
                                                                            else
                                                                              hj=(0,Inf)
                                                                            end
                                                                    else
                                                                      coefj=NTuple{3,Float64}((-βj*quanj-cj,-αj*quanj-bj,-quanj))
                                                                      allrealrootintervalnewtonregulafalsi(coefj,resppj,ppj)
                                                                      resTupj=(unsafe_load(resppj,1),unsafe_load(resppj,2))
                                                                      posSolPlusj=filter((x) -> x >0.0 , resTupj)
                                                                            if length(posSolPlusj)>0
                                                                              hj=(max(posSolPlusj...,ft-simt,firstguessH),Inf)
                                                                            else
                                                                              hj=(0,Inf)
                                                                            end
                                                                    end
                                                        end
                                                if hj[1]<=hi[1]<=hj[2] || hi[1]<=hj[1]<=hi[2] #overlap
                                        
                                                          h=min(hj[2],hi[2] )
                                                          if h==Inf && (hj[1]!=0.0 || hi[1]!=0.0) # except case both  0 sols  h1=(0,Inf)
                                                            h=max(hj[1],hi[1],ft-simt,firstguessH) #ft-simt in case they are both ver small...elaborate on his later
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
                                               end
                 =#

                else                   
                                     
                                      #find positive zeros f=+-Δ
                                    
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


                                          break
                                        else
                                          if resj[end][1]<resi[end][1]#remove last seg of resi
                                            resi=resi[1:end-1]
                                          else#remove last seg of resj
                                            resj=resj[1:end-1]
                                          end
                                        end

                                      end
                                      

               end

                                    
            
       
    
     
   
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
      
        tq[j]=simt 
      end #end second dependecy check
   end # end outer dependency check
   return iscycle
end   =#
