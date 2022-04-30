
function computeDerivative(index::Int, order::Int, jacobian::SVector{T,SVector{T,Float64}} ,x::MVector{O,Float64} , q::MVector{O,Float64}, tx::MVector{T,Float64}, tq::MVector{T,Float64}  )where{T,O}
  #for now i will not consider time ie derx=f(x)   
  x[(order+1)*index]=0
  for j = 1:T
   # x[(order+1)*index] += jacobian[j,index] * q[(order+1)*j-order] # use jacobian[j,index] when jacobian=transpose(jacobian)
    #x[(order+1)*index] += jacobian[index+(j-1)*T] * q[(order+1)*j-order] 
    x[(order+1)*index] += jacobian[index][j] * q[(order+1)*j-order] 
  end
end
function integrateState(j::Int,::Val{1},elapsed::Float64,x::MVector{O,Float64})where{O} 
      x[2*j-1]=x[2*j-1]+ elapsed * x[2*j]

  
end
function createDependencyMatrix(jacobian::SVector{N,SVector{N,Float64}}) where{N}
  # extract dependency matrix from jac
  #epselon=1e-6
#=   nRows=size(jac,1)
  nColumns=size(jac,2) =#
  #dep=SVector{nColumns,Array{Int}}([],[])###################################optimize meeeeeeee
  #dep=Vector{Array{Int}}(undef, N) 
#=  arr=[]
  for j=1:N
    tempvect=[]
       for i=1:N
           #if jac[i,j] < -epselon ||  jac[i,j] > epselon # different than zero
           #if jac[i+(j-1)*N]!=0
           if jacobian[i][j]!=0
               push!(tempvect,i)
           end
       end
       v=tuple(tempvect...)
       push!(arr,v) 
        
   end
   dep=(tuple(arr...)) =#
   #dep=Vector{Array{Int}}(undef, N)     ##########optimize ************
   pos = zeros(SVector{N,SVector{N,Int}})
   for j=1:N
        dep=[]
        for i=1:N
          if jacobian[i][j]!=0
                push!(dep,i)
          else
            push!(dep,0)
          end
        end
        pos=setindex(pos,dep,j)
    end
   
   
  return pos 
  #return dep 
  end
#=   epselon=1e-6
  nRows=size(jac,1)
  nColumns=size(jac,2)
  dep=Vector{Array{Int}}(undef, nColumns)     ##########optimize ************

  for j=1:nColumns
       dep[j]=Array{Float64}[]
       for i=1:nRows
           if jac[i,j] < -epselon ||  jac[i,j] > epselon # different than zero
               push!(dep[j],i)
           end
       end
   end
   
   #dep=SVector{nColumns,Array{Int}}([],[])

   return dep
end =#
function computeStates(initcond::SVector{N,Float64} ) where{N}
   #return length(initcond)
   return N
   #return 2

end
function getOrderfromSolverMethod(::Val{1})
    1
end
function getOrderfromSolverMethod(::Val{2})
  2
end
function getOrderfromSolverMethod(::Val{3})
  3
end
function modifyJacobian(jacobian :: SMatrix{N,N,Float64})where{N} #well remove smatrix from begining yay!!!!
  #test use of svector of svectos
  pos = zeros(SVector{N,SVector{N,Float64}})
  for i=1:N
    pos=setindex(pos,jacobian[i, :],i)
  end
  pos
end
#=
#without saving
function computeDerivative(states ::Int,index ::Int,order::Int,qssmodel :: QSSmodel,x ::  Vector{Float64} ,q::Vector{Float64} ,tx:: Vector{Float64} ,tq::Vector{Float64} )
   #for now i will not consider time ie derx=f(x)   
    der=0  
  # if length(x[(order+1)*index])!=0
   #       pop!(x[(order+1)*index])  
   #end
   for j = 1:states

       #lastQ=last(q[(order+1)*j-order])
       #der+=matrixEntry*lastQ
         der+=qssmodel.jacobian[index,j] * q[(order+1)*j-order]
   end
     x[(order+1)*index]=der   
end
=#

#=function computeDerivative(states::Int, index::Int, order::Int, jacobian::SMatrix{2,2,Float64} , x::MVector{4,Float64}, q::MVector{4,Float64}, tx::MVector{2,Float64}, tq::MVector{2,Float64})
  #for now i will not consider time ie derx=f(x)   
  #der = 0
x[(order+1)*index]=0
  for j = 1:states
    x[(order+1)*index] += jacobian[index, j] * q[(order+1)*j-order]
  end
  #x[(order+1)*index]= der


end=#
