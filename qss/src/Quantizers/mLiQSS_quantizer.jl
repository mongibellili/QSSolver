
#######################################################################################################################################################
function updateOtherApprox(::Val{1},j::Int,index::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=(q[index][0]-qaux[index][1])
   # @show q[index][0],qaux[index][1]
   # @show x[j][1],olddx[j][1]
   # println("aji before updateoher= ",a[j][index])
    if diffQ!=0
     a[j][index]=(x[j][1]-olddx[j][1])/diffQ
    else
     a[j][index]=0.0
    end
   # @show a[j][j],a[j][index]
   # println("u inside updateOther before update= ",u[j][index])
    u[j][index][1]=x[j][1]-a[j][j]*q[j][0]-a[j][index]*q[index][0]
   # println("u inside updateOther after update= ",u[j][index])
    return nothing
end
function updateOtherApprox(::Val{2},j::Int,index::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=q[index][0]-qaux[index][1]

  #  @show q[index][0],qaux[index][1]
  #  @show x[j][1],olddx[j][1]
  #  println("aji before updateoher= ",a[j][index])
    if diffQ != 0.0
        a[j][index]=(x[j][1]-olddx[j][1])/diffQ
    else
        a[j][index]=0.0
    end
  #  @show a[j][j],a[j][index]

  # println("u inside updateOther before update= ",u[j][index])
    u[j][index][1]=x[j][1]-a[j][j]*q[j][0]-a[j][index]*q[index][0]
    u[j][index][2]=2*x[j][2]-a[j][j]*q[j][1]-a[j][index]*q[index][1]
  #  println("u inside updateOther after update= ",u[j][index])
    #tu[index]=simt  # comment did nothing but it makes sense to keep it because more accurate since u is changed
    return nothing
end

#############################################################################################################################
function isCycle_and_simulUpdate(::Val{1},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tx::MVector{T,Float64},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64)where{T,O}
          elapsed = simt - tx[j]
          xjaux = x[j][0]+elapsed*x[j][1]
         dxj=a[j][index]*q[index][0]+a[j][j]*q[j][0]+u[j][index][1]
         iscycle=false
          if dxj*x[j][1]<0
            qjplus=xjaux+sign(dxj)*quantum[j]
            dxi=a[index][index]*q[index][0]+a[index][j]*qjplus+u[index][j][1]
            if dxi*x[index][1]<0
                iscycle=true
              qaux[j][1]=q[j][0]
              olddx[j][1]=x[j][1]
              aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];uij=u[index][j][1];uji=u[j][index][1];xi=x[index][0];xj=x[j][0];x1i=x[index][1];x1j=x[j][1]
              h = ft-simt
              Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
              qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
              qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
              if (abs(qi - xi) > 2*quantum[index] || abs(qj - xjaux) > 2*quantum[j]) # removing this did nothing...check @btime later
                h1 = (abs(quantum[index] / x1i));h2 = (abs(quantum[j] / x1j));
                h=min(h1,h2)
                Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
                if Δ==0
                  Δ=1e-12
                end
                qi = ((1-h*ajj)*(xi+h*uij)+h*aij*(xjaux+h*uji))/Δ
                qj = ((1-h*aii)*(xjaux+h*uji)+h*aji*(xi+h*uij))/Δ
              end
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
 println("nothing")
  xi1=x[index][1];xi2=2*x[index][2];xj1=x[j][1];xj2=2*x[j][2]
  e1 = simt - tu[j]
  xjaux = x[j](e1)# xAUX instead
  #println("ujj before update= ", u[j][j])
    u[j][j][1]=u[j][j][1]+e1*u[j][j][2]
   # println("ujj after update= ", u[j][j])
    e2 = simt - tq[j]
    tu[j]=simt
    tq[j]=simt
    #= qaux[j][1]=q[j][0] #used to be after if check
    olddx[j][1]=x[j][1] =#
    qaux[j][1]=q[j][0]+e2*q[j][1]
    olddx[j][1]=x[j][1]

    u[j][index][1]=u[j][j][1]-a[j][index]*qaux[index][1]
   dxj=a[j][index]*q[index][0]+a[j][j]*qaux[j][1]+u[j][index][1]
  # u[j][index][2]=u[j][j][2]-a[j][index]*qaux[index][2]#########this should be already updated... this formula brings nothing new
   ddxj=a[j][index]*q[index][1]+a[j][j]*q[j][1]+u[j][index][2]
   iscycle=false
  #=  @show a,u
   @show dxj,xj1
   @show ddxj,xj2 =#
    if (abs(dxj-x[j][1])>abs(dxj+x[j][1])/2 || abs(ddxj-2*x[j][2])>abs(ddxj+2*x[j][2])/2)
      #= if 11>simt>10
      @show j,dxj,x[j][1]
      end =#
      qjplus=xjaux-sign(ddxj)*quantum[j]
      h=sqrt(2*quantum[j]/abs(ddxj))
      dqjplus=(a[j][index]*(q[index][0]+h*q[index][1])+a[j][j]*qjplus+u[j][index][1]+h*u[j][index][2])/(1-h*a[j][j])
     # u[index][j][1]=u[index][index][1]-a[index][j]*qaux[j][1]
     # u[index][j][1]=u[index][index][1]-a[index][j]*q[j][0]
      dxi=a[index][index]*q[index][0]+a[index][j]*qjplus+u[index][j][1]
    #  u[index][j][2]=u[index][index][2]-a[index][j]*qaux[j][2]#########this should be already updated... this formula brings nothing new
      ddxi=a[index][index]*q[index][1]+a[index][j]*dqjplus+u[index][j][2]
     
      
      if (abs(dxi-x[index][1])>abs(dxi+x[index][1])/2 || abs(ddxi-2*x[index][2])>abs(ddxi+2*x[index][2])/2)
          iscycle=true
          @show dxi,xi1
          @show ddxi,xi2
        aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];uij=u[index][j][1];uji=u[j][index][1];uij2=u[index][j][2];uji2=u[j][index][2];xi=x[index][0];xj=x[j][0]
       # @show aii,ajj,aij,aji
       # @show uij,uji
        h = ft-simt
        Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
        α=(aii*(1-h*ajj)+h*aij*aji)/Δ
        β=(aij*(1-h*ajj)+h*aij*ajj)/Δ
        γ=(aji*(1-h*aii)+h*aji*aii)/Δ
        λ=(ajj*(1-h*aii)+h*aij*aji)/Δ
        α1=1+h*α-h*aii-h*h*(α*aii+γ*ajj)/2
        β1=h*(β-aij)-h*h*(β*aii+λ*aij)/2
        γ1=h*(γ-aji)-h*h*(α*aji+γ*ajj)/2
        λ1=1+h*(λ-ajj)-h*h*(β*aji+λ*ajj)/2
        θ10=(-h*(1-h*ajj)*(uij+h*uij2)-h*h*aij*(uji+h*uji2))/Δ
        θ11=xi+h*uij+h*h*(aii*(uij+h*uij2)+aij*(uji+h*uji2))/2
        θ20=(-h*(1-h*aii)*(uji+h*uji2)-h*h*aji*(uij+h*uij2))/Δ
        θ21=xjaux+h*uji+h*h*(ajj*(uji+h*uji2)+aji*(uij+h*uij2))/2
        θ12=θ10+θ11+uij2
        θ22=θ20+θ21+uji2
        Δ1=γ1*β1-α1*λ1
        qi=(λ1*θ12-β1*θ22)/Δ1
        qj=(-γ1*θ12-α1*θ22)/Δ1
        if (abs(qi - xi) > 2*quantum[index] || abs(qj - xjaux) > 2*quantum[j]) # removing this did nothing...check @btime later
          h1 = sqrt(abs(quantum[index]/x[index][2]));h2 = sqrt(abs(quantum[j]/x[j][2]));
          h=min(h1,h2)
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          α=(aii*(1-h*ajj)+h*aij*aji)/Δ
          β=(aij*(1-h*ajj)+h*aij*ajj)/Δ
          γ=(aji*(1-h*aii)+h*aji*aii)/Δ
          λ=(ajj*(1-h*aii)+h*aij*aji)/Δ
          α1=1+h*α-h*aii-h*h*(α*aii+γ*ajj)/2
          β1=h*(β-aij)-h*h*(β*aii+λ*aij)/2
          γ1=h*(γ-aji)-h*h*(α*aji+γ*ajj)/2
          λ1=1+h*(λ-ajj)-h*h*(β*aji+λ*ajj)/2
          θ10=(-h*(1-h*ajj)*(uij+h*uij2)-h*h*aij*(uji+h*uji2))/Δ
          θ11=xi+h*uij+h*h*(aii*(uij+h*uij2)+aij*(uji+h*uji2))/2
          θ20=(-h*(1-h*aii)*(uji+h*uji2)-h*h*aji*(uij+h*uij2))/Δ
          θ21=xjaux+h*uji+h*h*(ajj*(uji+h*uji2)+aji*(uij+h*uij2))/2
          θ12=θ10+θ11+uij2
          θ22=θ20+θ21+uji2
          Δ1=γ1*β1-α1*λ1
          qi=(λ1*θ12-β1*θ22)/Δ1
          qj=(-γ1*θ12-α1*θ22)/Δ1
        end
        maxIter=1000
        while (abs(qi - xi) > 2*quantum[index] || abs(qj - xjaux) > 2*quantum[j]) && (maxIter>0)
          maxIter-=1
          h1 = h * sqrt(2*quantum[index] / abs(qi - xi));
          h2 = h * sqrt(2*quantum[j] / abs(qj - xjaux));
          h=min(h1,h2)
          Δ=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
          α=(aii*(1-h*ajj)+h*aij*aji)/Δ
          β=(aij*(1-h*ajj)+h*aij*ajj)/Δ
          γ=(aji*(1-h*aii)+h*aji*aii)/Δ
          λ=(ajj*(1-h*aii)+h*aij*aji)/Δ
          α1=1+h*α-h*aii-h*h*(α*aii+γ*ajj)/2
          β1=h*(β-aij)-h*h*(β*aii+λ*aij)/2
          γ1=h*(γ-aji)-h*h*(α*aji+γ*ajj)/2
          λ1=1+h*(λ-ajj)-h*h*(β*aji+λ*ajj)/2
          θ10=((1-h*ajj)*(uij+h*uij2)+h*aij*(uji+h*uji2))/Δ
          θ11=xi+h*uij+h*h*(aii*(uij+h*uij2)+aij*(uji+h*uji2))/2
          θ20=((1-h*aii)*(uji+h*uji2)+h*aji*(uij+h*uij2))/Δ
          θ21=xjaux+h*uji+h*h*(ajj*(uji+h*uji2)+aji*(uij+h*uij2))/2
          θ12=-h*θ10+θ11+uij2
          θ22=-h*θ20+θ21+uji2
          Δ1=γ1*β1-α1*λ1
          qi=(λ1*θ12-β1*θ22)/Δ1
          qj=(-γ1*θ12-α1*θ22)/Δ1
        end
        @show maxIter
       #println("after while loop")
        #@show qi,xi,quantum[index]
       # @show qj,xjaux,quantum[j]
       # limitedPrint=2
        q[index][0]=qi# store back helper vars
        q[j][0]=qj
        q[index][1]=α*qi+β*qj+θ10
        q[j][1]=γ*qi+λ*qj+θ20
      end #end second dependecy check
  end # end outer dependency check
  return iscycle
end