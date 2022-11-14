function isCycle_and_simulUpdate(::Val{2},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
    aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];xi=x[index][0];xj=x[j][0];qi=q[index][0];qj=q[j][0];qi1=q[index][1];qj1=q[j][1]
    quanj=quantum[j]
    quani=quantum[index]
    xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
    e1 = simt - tx[j]
    xjaux = x[j](e1)# xAUX instead
     #e3=simt - tu[j]
    u[j][j][1]=u[j][j][1]+e1*u[j][j][2]  #e3 does not work here
    e2 = simt - tq[j]
     # tx[j]=simt   # do not update tx[j] here ie do not uncomment because we are not doing a real update it is just a prediction
     # tu[j]=simt  # does not matter
     # tq[j]=simt   # does not matter
    qaux[j][1]=qj+e2*qj1  ###################should i eliminate the time update?????????
    olddx[j][1]=xj1
  
    u[j][index][1]=u[j][j][1]-aji*qaux[index][1]# using q[i][0] creates a really huge bump at 18 (no go)
    uji=u[j][index][1]
    dxj=aji*qi+ajj*qaux[j][1]+uji
     #  u[j][index][2]=u[j][j][2]-a[j][index]*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
    u[j][index][2]=u[j][j][2]-ajj*qaux[index][1] # more cycles ...shaky with no bumps
    uji2=u[j][index][2]
    ddxj=aji*qi1+ajj*qj1+uji2
    iscycle=false
    if (abs(dxj-xj1)>(abs(dxj+xj1)/2) || abs(ddxj-2*xj2)>(abs(ddxj+2*xj2)/2))    
      qjplus=xjaux-sign(ddxj)*quanj
      h=sqrt(quanj/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
      dqjplus=(aji*(qi+h*qi1)+ajj*qjplus+uji+h*uji2)/(1-h*ajj)
      u[index][j][1]=u[index][index][1]-aij*qaux[j][1]
      #u[index][j][1]=u[index][index][1]-a[index][j]*q[j][0]  # shifts down at 18
      uij=u[index][j][1]
      dxi=aii*qi+aij*qjplus+uij
      u[index][j][2]=u[index][index][2]-aij*qj1#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
      uij2=u[index][j][2]
      ddxi=aii*qi1+aij*dqjplus+uij2
      if (abs(dxi-xi1)>(abs(dxi+xi1)/2) || abs(ddxi-xi2)>(abs(ddxi+xi2)/2))
          iscycle=true
          A=[aii aij;aji ajj]
          I=[1 0;0 1]
          U=[uij;uji]
          U2=[uij2;uji2]
          X=[xi;xjaux]
          h = ft-simt
          N=inv(I-h*A)
          Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2)        
          qi=Q[1]
          qj=Q[2]
          if (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) # removing this did nothing...check @btime later
            h1 = sqrt(abs(quani/xi2));h2 = sqrt(abs(quanj/xj2));
            h=min(h1,h2)
            N=inv(I-h*A)
            Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2)         
            qi=Q[1]
            qj=Q[2]
          end
          maxIter=10000
          while (abs(qi - xi) > 2*quani || abs(qj - xjaux) > 2*quanj) && (maxIter>0)
            maxIter-=1
            h1 = h * sqrt(2*quani / abs(qi - xi));
            h2 = h * sqrt(2*quanj / abs(qj - xjaux));
            h=min(h1,h2)
            N=inv(I-h*A)
            Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2)         
            qi=Q[1]
            qj=Q[2]
          end
  
          q[index][0]=qi# store back helper vars
          q[j][0]=qj
          Q1=N*(A*Q+U+h*U2)
          q[index][1]=Q1[1]# store back helper vars
          q[j][1]=Q1[2]
   #=  h = ft-simt
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          α=(aii*(1-h*ajj)+h*aij*aji)/Δ
          β=(aij*(1-h*ajj)+h*aij*ajj)/Δ
          γ=(aji*(1-h*aii)+h*aji*aii)/Δ
          λ=(ajj*(1-h*aii)+h*aij*aji)/Δ
          α1=1+h*α-h*aii-h*h*(α*aii+γ*ajj)/2
          β1=h*(β-aij)-h*h*(β*aii+λ*aij)/2
          γ1=h*(γ-aji)-h*h*(α*aji+γ*ajj)/2
          λ1=1+h*(λ-ajj)-h*h*(β*aji+λ*ajj)/2
          z1=aii*(1-h*ajj)+h*aji*aij
          z2=aii*h*aij+aij*(1-h*aii)
          z3=aji*(1-h*ajj)+h*ajj*aji
          z4=aji*h*aij+ajj*(1-h*aii)
          θ10=(-h*(1-h*ajj)*(uij+h*uij2)-h*h*aij*(uji+h*uji2))/Δ
          θ11=xi+h*uij+h*h*(z1*(uij+h*uij2)+z2*(uji+h*uji2))/2
          θ20=(-h*(1-h*aii)*(uji+h*uji2)-h*h*aji*(uij+h*uij2))/Δ
          θ21=xjaux+h*uji+h*h*(z3*(uji+h*uji2)+z4*(uij+h*uij2))/2
          θ12=θ10+θ11+uij2
          θ22=θ20+θ21+uji2
          Δ1=γ1*β1-α1*λ1
          qi=(λ1*θ12-β1*θ22)/Δ1
          qj=(-γ1*θ12-α1*θ22)/Δ1
          println("lets see qi and qj after initial change using ft: ")
          @show qi,qj
          
          
          if (abs(qi - xi) > 2*quantum[index] || abs(qj - xjaux) > 2*quantum[j]) # removing this did nothing...check @btime later
            h1 = sqrt(abs(quantum[index]/x[index][2]));h2 = sqrt(abs(quantum[j]/x[j][2]));
            @show h1,h2
            h=min(h1,h2)
            println("h inside sqrt(quan/ddx)= ",h)
  
            Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
            α=(aii*(1-h*ajj)+h*aij*aji)/Δ
            β=(aij*(1-h*ajj)+h*aij*ajj)/Δ
            γ=(aji*(1-h*aii)+h*aji*aii)/Δ
            λ=(ajj*(1-h*aii)+h*aij*aji)/Δ
            α1=1+h*α-h*aii-h*h*(α*aii+γ*ajj)/2
            β1=h*(β-aij)-h*h*(β*aii+λ*aij)/2
            γ1=h*(γ-aji)-h*h*(α*aji+γ*ajj)/2
            λ1=1+h*(λ-ajj)-h*h*(β*aji+λ*ajj)/2
            z1=aii*(1-h*ajj)+h*aji*aij
            z2=aii*h*aij+aij*(1-h*aii)
            z3=aji*(1-h*ajj)+h*ajj*aji
            z4=aji*h*aij+ajj*(1-h*aii)
            θ10=(-h*(1-h*ajj)*(uij+h*uij2)-h*h*aij*(uji+h*uji2))/Δ
            θ11=xi+h*uij+h*h*(z1*(uij+h*uij2)+z2*(uji+h*uji2))/2
            θ20=(-h*(1-h*aii)*(uji+h*uji2)-h*h*aji*(uij+h*uij2))/Δ
            θ21=xjaux+h*uji+h*h*(z3*(uji+h*uji2)+z4*(uij+h*uij2))/2
            θ12=θ10+θ11+uij2
            θ22=θ20+θ21+uji2
            Δ1=γ1*β1-α1*λ1
            qi=(λ1*θ12-β1*θ22)/Δ1
            qj=(-γ1*θ12-α1*θ22)/Δ1
          end
          println("lets see qi and qj after second change using sqrt(quan/x2): ")
          @show qi,qj
        
          maxIter=10000
          while (abs(qi - xi) > 2*quantum[index] || abs(qj - xjaux) > 2*quantum[j]) && (maxIter>0)
            maxIter-=1
            h1 = h * sqrt(2*quantum[index] / abs(qi - xi));
            h2 = h * sqrt(2*quantum[j] / abs(qj - xjaux));
            h=min(h1,h2)
            Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
            if Δ==0
              Δ=1e-12
            end
            α=(aii*(1-h*ajj)+h*aij*aji)/Δ
            β=(aij*(1-h*ajj)+h*aij*ajj)/Δ
            γ=(aji*(1-h*aii)+h*aji*aii)/Δ
            λ=(ajj*(1-h*aii)+h*aij*aji)/Δ
            α1=1+h*α-h*aii-h*h*(α*aii+γ*ajj)/2
            β1=h*(β-aij)-h*h*(β*aii+λ*aij)/2
            γ1=h*(γ-aji)-h*h*(α*aji+γ*ajj)/2
            λ1=1+h*(λ-ajj)-h*h*(β*aji+λ*ajj)/2
            z1=aii*(1-h*ajj)+h*aji*aij
            z2=aii*h*aij+aij*(1-h*aii)
            z3=aji*(1-h*ajj)+h*ajj*aji
            z4=aji*h*aij+ajj*(1-h*aii)
            θ10=(-h*(1-h*ajj)*(uij+h*uij2)-h*h*aij*(uji+h*uji2))/Δ
            θ11=xi+h*uij+h*h*(z1*(uij+h*uij2)+z2*(uji+h*uji2))/2
            θ20=(-h*(1-h*aii)*(uji+h*uji2)-h*h*aji*(uij+h*uij2))/Δ
            θ21=xjaux+h*uji+h*h*(z3*(uji+h*uji2)+z4*(uij+h*uij2))/2
            θ12=-h*θ10+θ11+uij2
            θ22=-h*θ20+θ21+uji2
            Δ1=γ1*β1-α1*λ1
            if Δ1==0
              Δ1=1e-12
            end
            qi=(λ1*θ12-β1*θ22)/Δ1
            qj=(-γ1*θ12-α1*θ22)/Δ1
          end
          @show maxIter
         println("after while loop")
          @show qi,xi,quantum[index]
          @show qj,xjaux,quantum[j]
         # limitedPrint=2
          q[index][0]=qi# store back helper vars
          q[j][0]=qj
          q[index][1]=α*qi+β*qj+θ10
          q[j][1]=γ*qi+λ*qj+θ20 
   =#
        end #end second dependecy check
    end # end outer dependency check
    return iscycle
  end