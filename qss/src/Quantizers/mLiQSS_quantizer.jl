
#######################################################################################################################################################
function updateOtherApprox(::Val{1},j::Int,index::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=(q[index][0]-qaux[index][1])
   # @show q[index][0],qaux[index][1]
   # @show x[j][1],olddx[j][1]
   # println("aji before updateoher= ",a[j][index])
   #if 1.21<simt<1.226597
   #@show simt 
   #@show j
   # println("aji before updateoher= ",a[j][index])
   # @show q[index][0],qaux[index][1]
   # @show x[j][1],olddx[j][1]
  # end
    if diffQ!=0
     a[j][index]=(x[j][1]-olddx[j][1])/diffQ
    else
     a[j][index]=0.0
    end
   # if 1.21<simt<1.226597
   #   println("aji afterupdateoher= ",a[j][index])
      
    # end
   # @show a[j][j],a[j][index]
   # println("u inside updateOther before update= ",u[j][index])
    u[j][index][1]=x[j][1]-a[j][j]*q[j][0]-a[j][index]*q[index][0]
   # println("u inside updateOther after update= ",u[j][index])
    return nothing
end
function updateOtherApprox(::Val{2},j::Int,index::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=q[index][0]-qaux[index][1]
   #if 1.54<simt<1.546597
   
    
    
  # end
   # println("aji before updateoher= ",a[j][index])
    if diffQ != 0.0
        a[j][index]=(x[j][1]-olddx[j][1])/diffQ
    else
        a[j][index]=0.0
    end
  #  @show a[j][j],a[j][index]
 # if 1.54<simt<1.546597
   # println("a$j$i afterupdateoher= ",a[j][index])
    
  # end
  # println("u inside updateOther before update= ",u[j][index])
    u[j][index][1]=x[j][1]-a[j][j]*q[j][0]-a[j][index]*q[index][0]
    u[j][index][2]=2*x[j][2]-a[j][j]*q[j][1]-a[j][index]*q[index][1]
  #  println("u inside updateOther after update= ",u[j][index])
    #tu[index]=simt  # comment did nothing but it makes sense to keep it because more accurate since u is changed
    return nothing
end

#############################################################################################################################
function isCycle_and_simulUpdate(::Val{1},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
  aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];uij=u[index][j][1];uji=u[j][index][1];xi=x[index][0];xj=x[j][0];x1i=x[index][1];x1j=x[j][1]
  qi=q[index][0];qj=q[j][0]
  
  qaux[j][1]=q[j][0]
              olddx[j][1]=x[j][1]
        elapsed = simt - tx[j]
          xjaux = x[j][0]+elapsed*x[j][1]
         dxj=a[j][index]*q[index][0]+a[j][j]*q[j][0]+u[j][index][1]
         iscycle=false

         #= if 1.53<simt < 1.5429667149224133
          @show simt
          @show aii,ajj,aij,aji
      @show uij,uji
      @show qi,qj
      @show xi,xjaux
        end =#

          if dxj*x[j][1]<0
            qjplus=xjaux+sign(dxj)*quantum[j]
            dxi=a[index][index]*q[index][0]+a[index][j]*qjplus+u[index][j][1]
            if dxi*x[index][1]<0
                iscycle=true
              #= qaux[j][1]=q[j][0]
              olddx[j][1]=x[j][1] =#
              
       @show simt
              @show aii,ajj,aij,aji
        @show uij,uji
        @show qi,qj
        @show xi,xjaux
              
              h = ft-simt
              Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
              qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
              qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
              if (abs(qi - xi) > 2*quantum[index] || abs(qj - xjaux) > 2*quantum[j]) # removing this did nothing...check @btime later
                h1 = (abs(quantum[index] / x1i));h2 = (abs(quantum[j] / x1j));
                h=min(h1,h2)
                @show h1,h2
               
                println("h inside sqrt(quan/ddx)= ",h)
                Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                if Δ==0
                  Δ=1e-12
                end
                qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
                qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
              end
              println("lets see qi and qj after second change using sqrt(quan/x2): ")
              @show qi,qj
              maxIter=100
              while (abs(qi - xi) > 2*quantum[index] || abs(qj - xjaux) > 2*quantum[j]) && (maxIter>0)
                maxIter-=1
                h1 = h * (2*quantum[index] / abs(qi - xi));
                h2 = h * (2*quantum[j] / abs(qj - xjaux));
                h=min(h1,h2)
                Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                if Δ==0
                 # println("delta==0")
                  Δ=1e-12
                end
                qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
                qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
              end
              @show maxIter
              println("after while loop")
               @show qi,xi,quantum[index]
               @show qj,xjaux,quantum[j]
             #println("after while loop")
              #@show qi,xi,quantum[index]
             # @show qj,xjaux,quantum[j]
             # limitedPrint=2
              q[index][0]=qi# store back helper vars
              q[j][0]=qj
            end #end second dependecy check
        end # end outer dependency check
        return iscycle
end
         

function isCycle_and_simulUpdate(::Val{2},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
 
  xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
# println("nothing")

  xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  e1 = simt - tx[j]
  xjaux = x[j](e1)# xAUX instead
 # println("ujj before update= ", u[j][j])
 e3=simt - tu[j]
    u[j][j][1]=u[j][j][1]+e1*u[j][j][2]  #e3 does not work here
   # println("ujj after update= ", u[j][j])
    e2 = simt - tq[j]
   # tx[j]=simt   # do not update tx[j] here ie do not uncomment because we are not doing a real update it is just a prediction
   # tu[j]=simt  # does not matter
   # tq[j]=simt   # does not matter
    #= qaux[j][1]=q[j][0] #used to be after if check
    olddx[j][1]=x[j][1] =#
    qaux[j][1]=q[j][0]+e2*q[j][1]  ###################should i eliminate the time update?????????
    olddx[j][1]=x[j][1]
   ujitemp=u[j][index][1]
   uji2temp=u[j][index][2]
   u[j][index][1]=u[j][j][1]-a[j][index]*qaux[index][1]# using q[i][0] creates a really huge bump at 18 (no go)
   
   dxj=a[j][index]*q[index][0]+a[j][j]*qaux[j][1]+u[j][index][1]
 #  u[j][index][2]=u[j][j][2]-a[j][index]*qaux[index][2]#less cycles but with a bump at 1.5...ft20: smooth with some bumps
   u[j][index][2]=u[j][j][2]-a[j][j]*qaux[index][1] # more cycles ...shaky with no bumps
   ddxj=a[j][index]*q[index][1]+a[j][j]*q[j][1]+u[j][index][2]
   iscycle=false

  #=  @show a,u
   @show dxj,xj1
   @show ddxj,xj2 =#
    if (abs(dxj-x[j][1])>(abs(dxj+x[j][1])/2) || abs(ddxj-2*x[j][2])>(abs(ddxj+2*x[j][2])/2))
      #= if 11>simt>10
      @show j,dxj,x[j][1]
      end =#
      
      qjplus=xjaux-sign(ddxj)*quantum[j]
      h=sqrt(quantum[j]/abs(ddxj))#2*quantum funny oscillating graph; xj2 vibrating
      dqjplus=(a[j][index]*(q[index][0]+h*q[index][1])+a[j][j]*qjplus+u[j][index][1]+h*u[j][index][2])/(1-h*a[j][j])
      u[index][j][1]=u[index][index][1]-a[index][j]*qaux[j][1]
      #u[index][j][1]=u[index][index][1]-a[index][j]*q[j][0]  # shifts down at 18
      dxi=a[index][index]*q[index][0]+a[index][j]*qjplus+u[index][j][1]
     u[index][j][2]=u[index][index][2]-a[index][j]*q[j][1]#########qaux[j][2] updated in normal Qupdate..ft=20 slightly shifts up
      ddxi=a[index][index]*q[index][1]+a[index][j]*dqjplus+u[index][j][2]
     
      
      if (abs(dxi-x[index][1])>(abs(dxi+x[index][1])/2) || abs(ddxi-2*x[index][2])>(abs(ddxi+2*x[index][2])/2))
          iscycle=true
          
        #  @show q[index][0],q[j][0]
       #  @show x[index][0],xjaux
       #   @show abs(q[index][0] - x[index][0]) > 2*quantum[index] ,2*quantum[index] 
        #  @show abs(q[j][0] - xjaux) > 2*quantum[j],2*quantum[j]
        # @show j
         #=  @show simt
          @show dxi,xi1
          @show ddxi,xi2
          @show dxj,xj1
          @show ddxj,xj2 =#


        aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];uij=u[index][j][1];uji=u[j][index][1];uij2=u[index][j][2];uji2=u[j][index][2];xi=x[index][0];xj=x[j][0]
        qi=q[index][0];qj=q[j][0]
        #ajj=-1;aii=-1;aji=87;aij=3.01;uji=2.61;uij=3.68
        #= @show simt
        @show aii,ajj,aij,aji
        @show uij,uji
        @show qi,qj
        @show xi,xjaux =#
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
        if (abs(qi - xi) > 2*quantum[index] || abs(qj - xjaux) > 2*quantum[j]) # removing this did nothing...check @btime later
          h1 = sqrt(abs(quantum[index]/x[index][2]));h2 = sqrt(abs(quantum[j]/x[j][2]));
        #  @show h1,h2
          h=min(h1,h2)
          N=inv(I-h*A)
          Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2)
          
          qi=Q[1]
          qj=Q[2]
        end
        maxIter=10000
        while (abs(qi - xi) > 2*quantum[index] || abs(qj - xjaux) > 2*quantum[j]) && (maxIter>0)
          maxIter-=1
          h1 = h * sqrt(2*quantum[index] / abs(qi - xi));
          h2 = h * sqrt(2*quantum[j] / abs(qj - xjaux));
          h=min(h1,h2)
          N=inv(I-h*A)
          Q=inv(I-h*A+h*N*A-h*h*A*N*A/2)*(((h*h/2)*A-h*I)*N*(U+h*U2)+X+h*U+h*h*U2)
          
          qi=Q[1]
          qj=Q[2]
        end
       # @show maxIter
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
       #   @show a
       #   @show u
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
        q[j][1]=γ*qi+λ*qj+θ20 =#

        
      end #end second dependecy check
  end # end outer dependency check
  return iscycle
end