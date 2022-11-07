
#######################################################################################################################################################
function updateOtherApprox(::Val{1},j::Int,index::Int,x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tu::MVector{T,Float64},simt::Float64)where{T,O}
    diffQ=(q[index][0]-qaux[index][1])
    if diffQ!=0
     a[j][index]=(x[j][1]-olddx[j][1])/diffQ
    #=  if x[j][1]==olddx[j][1]
       println("x[$j][1]-olddx[$j][1]")
       @show olddx[j][1]
       @show qaux[j][1],qtemp
       @show qaux[index][1],q[index][0]
     end =#
    else
     a[j][index]=0.0
     #println("a[$j][$index]==0")
    end
    u[j][index][1]=x[j][1]-a[j][j]*q[j][0]-a[j][index]*q[index][0]
    return nothing
end
#############################################################################################################################
function isCycle_and_simulUpdate(::Val{1},index::Int,j::Int, x::Vector{Taylor0{Float64}},q::Vector{Taylor0{Float64}}, quantum::Vector{Float64},a::MVector{T,MVector{T,Float64}},u::MVector{T,MVector{T,MVector{O,Float64}}},qaux::MVector{T,MVector{O,Float64}},olddx::MVector{T,MVector{O,Float64}},tq::MVector{T,Float64},tu::MVector{T,Float64},simt::Float64,ft::Float64,elapsed::Float64)where{T,O}
    
          #xjaux = x[j](elapsed)# xAUX instead
          xjaux = x[j][0]+elapsed*x[j][1]
         # i=index
         dxj=a[j][index]*q[index][0]+a[j][j]*q[j][0]+u[j][index][1]
        # dxi=a[index][index]*qaux[index][1]+a[index][j]*q[j][0]+u[index][j][1]
        #=  if printcount%50==0
         @show a
         @show u
        # @show dxi,x[index][1]
         #@show q[index][0],q[j][0]
         @show dxj,x[j][1]
         end =#
         # u[j][index][1]=u[j][j][1]-a[j][index]*qaux[index][1]
         iscycle=false
          if dxj*x[j][1]<0
            #= if 11>simt>10
            @show j,dxj,x[j][1]
            end =#
            qjplus=xjaux+sign(dxj)*quantum[j]
           # u[index][j][1]=u[index][index][1]-a[index][j]*q[j][0]
            dxi=a[index][index]*q[index][0]+a[index][j]*qjplus+u[index][j][1]
            #= if printcount%50==0
               @show dxi,x[index][1]
            end =#
         #   
           #=  if 11>simt>10
              @show a[index][index],q[index][0],a[index][j],qjplus,u[index][j][1]
              @show index,dxi,x[index][1]
            end =#
            if dxi*x[index][1]<0
                iscycle=true
              qaux[j][1]=q[j][0]
              olddx[j][1]=x[j][1]
              aii=a[index][index];ajj=a[j][j];aij=a[index][j];aji=a[j][index];uij=u[index][j][1];uji=u[j][index][1];xi=x[index][0];xj=x[j][0];x1i=x[index][1];x1j=x[j][1]
             # @show aii,ajj,aij,aji
             # @show uij,uji
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
             