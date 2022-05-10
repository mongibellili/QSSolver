
function modifyJacobian(jacobian :: SMatrix{N,N,Float64})where{N} #well test removing smatrix from begining 
  pos = zeros(SVector{N,SVector{N,Float64}})
  for i=1:N
    pos=setindex(pos,jacobian[i, :],i)
  end
  pos
end

function createDependencyMatrix(jacobian::SVector{N,SVector{N,Float64}}) where{N}
   dep = zeros(SVector{N,SVector{N,Int}})
   for j=1:N
        arr=[]
        for i=1:N
          if jacobian[i][j]!=0#if jac[i,j] < -epselon ||  jac[i,j] > epselon # different than zero
                push!(arr,i)
          else
            push!(arr,0) # for some reason inner vectors want to have equal sizes... same when used tuples of tuples!!!
          end
        end
        dep=setindex(dep,arr,j)
    end  
    return dep  #dep=(tuple(arr...))  could use tuples; they also want equal size innner tuples!!!!!
  end

function computeStates(initcond::SVector{N,Float64} ) where{N}
   return N#return length(initcond)   
end
function computeInputs(u::SVector{S,Float64} ) where{S}
  return S#return length(initcond)   
end
function computeInputs(u::SVector{S,Function} ) where{S}
  return S#return length(initcond)   
end
function derivateJacobian(jacobian :: SMatrix{N,N,Float64})where{N}
return jacobian*jacobian #this is for  time-independent jacobians# run once: not needed cuz i used derQ instead
end

getOrderfromSolverMethod(::Val{1})=1
getOrderfromSolverMethod(::Val{2})=2
getOrderfromSolverMethod(::Val{3})=3
#getOrderfromSolverMethod(::Val{4})=1 #liqss1 #having more than 3 methods is causing an error 
#getOrderfromSolverMethod(::Val{5})=2  #liqss2