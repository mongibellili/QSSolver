using BenchmarkTools
function sort3(tup)
  h1=min(tup[1],tup[2],tup[3])
  h3=max(tup[1],tup[2],tup[3])
  if tup[1]!=h1 && tup[1]!=h3
    h2=tup[1]
  elseif tup[2]!=h1 && tup[2]!=h3
    h2=tup[2]
  else
    h2=tup[3]
  end
  res=(h1,h2,h3)
end
#= function constructIntrval1(tup)
  lenT=length(tup)
  #res=()
  if lenT==0
    res=((0,Inf),)
  elseif lenT==1
    res=((0,tup[1]),)
  elseif lenT==2
    h1=min(tup[1],tup[2])
    h4=max(tup[1],tup[2])
    res=((0,h1),(h4,Inf))
  elseif lenT==3
    h1,h2,h3=sort3(tup)
    res=((0,h1),(h2,h3))
  else#4
    h1=min(tup[1],tup[2],tup[3],tup[4])
    if tup[1]==h1
      h2,h3,h4=sort3((tup[2],tup[3],tup[4]))
    elseif tup[2]==h1
      h2,h3,h4=sort3((tup[1],tup[3],tup[4]))
    elseif tup[3]==h1
      h2,h3,h4=sort3((tup[2],tup[1],tup[4]))
    else
      h2,h3,h4=sort3((tup[2],tup[3],tup[1]))
    end
    res=((0,h1),(h2,h3),(h4,Inf))
  end

  return res
end =#
constructIntrval(tup::Tuple{})=((0,Inf),)
constructIntrval(tup::Tuple{Float64})=((0,tup[1]),)
function constructIntrval(tup::Tuple{Float64,Float64})
    h1=min(tup[1],tup[2])
    h4=max(tup[1],tup[2])
    res=((0,h1),(h4,Inf))
end
function constructIntrval(tup::NTuple{3, Float64})
    h1,h2,h3=sort3(tup)
    res=((0,h1),(h2,h3))
end
function constructIntrval(tup::NTuple{4, Float64})
    h1=min(tup[1],tup[2],tup[3],tup[4])
    if tup[1]==h1
      h2,h3,h4=sort3((tup[2],tup[3],tup[4]))
    elseif tup[2]==h1
      h2,h3,h4=sort3((tup[1],tup[3],tup[4]))
    elseif tup[3]==h1
      h2,h3,h4=sort3((tup[2],tup[1],tup[4]))
    else
      h2,h3,h4=sort3((tup[2],tup[3],tup[1]))
    end
    res=((0,h1),(h2,h3),(h4,Inf))
 
end


#hl=(tu...,tb...)
#@show hl
#= newtu=hl[1:end-1]
@show newtu =#
#= @show constructIntrval1(hl)
@show constructIntrval(hl)


@btime constructIntrval1(hl)
@btime constructIntrval(hl) =#

#@show res
#= tu=()
display(typeof(tu)) =#

tu=(2,)
@show tu[1:end-1]


function nmisCycle_and_simulUpdate(gh,#= simuldeltaiVals,simuldeltajVals,simulqxiVals,simulqxjVals, simulHTimes,simulHVals, =#::Val{1},index::Int,j::Int,dirI::Float64,firstguessH::Float64, x::Vector{Taylor0},q::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64)
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
                                      pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
                                    
                                      respp = pointer(Vector{Float64}(undef, 2))
                                      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2) #force negative values
                                      
                                      #find positive zeros f=+-Δ
                                      bi=aii*xi+aij*xj+uij;ci=aij*(aji*xi+uji)-ajj*(aii*xi+uij);αi=-ajj-aii;βi=aii*ajj-aij*aji
                                      #= coefi=@SVector [quani,αi*quani-bi,βi*quani-ci]#
                                      posSolPlusi= PosRoot(coefi, Val(2)) =#
                                      coefi=NTuple{3,Float64}((βi*quani-ci,αi*quani-bi,quani))
                                      allrealrootintervalnewtonregulafalsi(coefi,respp,pp)
                                      resTupi=(unsafe_load(respp,1),unsafe_load(respp,2))
                                      posSolPlusi=filter((x) -> x >0.0 , resTupi)
                                    #=  coefi=@SVector [-quani,-αi*quani-bi,-βi*quani-ci]
                                      posSolminusi= PosRoot(coefi, Val(2)) =#
                                      coefi2=NTuple{3,Float64}((-βi*quani-ci,-αi*quani-bi,-quani))
                                      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2) #force negative values
                                      allrealrootintervalnewtonregulafalsi(coefi2,respp,pp)
                                      resTupi2=(unsafe_load(respp,1),unsafe_load(respp,2))
                                      posSolminusi=filter((x) -> x >0.0 , resTupi2)
                                      
                                      posSoli=(posSolPlusi...,posSolminusi...)

                                      bj=ajj*xj+aji*xi+uji;cj=aji*(aij*xj+uij)-aii*(ajj*xj+uji);αj=-aii-ajj;βj=ajj*aii-aji*aij

                                      #= coefj=@SVector [quanj,αj*quanj-bj,βj*quanj-cj]#
                                      posSolPlusj= PosRoot(coefj, Val(2)) =#
                                      coefj=NTuple{3,Float64}((βj*quanj-cj,αj*quanj-bj,quanj))
                                      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2) #force negative values
                                      allrealrootintervalnewtonregulafalsi(coefj,respp,pp)
                                      resTupj=(unsafe_load(respp,1),unsafe_load(respp,2))
                                      posSolPlusj=filter((x) -> x >0.0 , resTupj)
                                      #= coefj2=@SVector [-quanj,-αj*quanj-bj,-βj*quanj-cj]#
                                      posSolminusj= PosRoot(coefj2, Val(2)) =#
                                      coefj2=NTuple{3,Float64}((-βj*quanj-cj,-αj*quanj-bj,-quanj))
                                      unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2) #force negative values
                                      allrealrootintervalnewtonregulafalsi(coefj2,respp,pp)
                                      resTupj2=(unsafe_load(respp,1),unsafe_load(respp,2))
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
end  