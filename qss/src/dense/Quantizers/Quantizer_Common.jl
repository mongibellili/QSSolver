@inline function integrateState(::Val{0}, x::Taylor0,elapsed::Float64) 
  #nothing: created for elapse-updating q in order1 which does not happen
end

@inline function integrateState(::Val{1}, x::Taylor0,elapsed::Float64) 
  x.coeffs[1] = x(elapsed)
end
@inline function integrateState(::Val{2}, x::Taylor0,elapsed::Float64) 
  x.coeffs[1] = x(elapsed)
  x.coeffs[2] = x.coeffs[2]+elapsed*x.coeffs[3]*2
end

#= @inline function integrateState(::Val{3}, x::Taylor0,elapsed::Float64) 
  x.coeffs[1] = x(elapsed)
  x.coeffs[2] = x.coeffs[2]+elapsed*x.coeffs[3]*2+elapsed*elapsed*x.coeffs[4]*3
  x.coeffs[3] = x.coeffs[3]+elapsed*x.coeffs[4]*3#(x.coeffs[3]*2+elapsed*x.coeffs[4]*6)/2
end =#


######################################################################################################################################"
function computeDerivative( ::Val{1}  ,x::Taylor0,f::Taylor0 )
    x.coeffs[2] =f.coeffs[1]
    return nothing
end
function computeDerivative( ::Val{2}  ,x::Taylor0,f::Taylor0) 
    x.coeffs[2] =f.coeffs[1]
    x.coeffs[3] =f.coeffs[2]/2
    return nothing
end

#= function computeDerivative( ::Val{3}  ,x::Taylor0,f::Taylor0)
  x.coeffs[2] =f[0]
  x.coeffs[3]=f.coeffs[2]/2
  x.coeffs[4]=f.coeffs[3]/3 # coeff3*2/6
  return nothing
end =#



function updateLinearApprox(i::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},simt)where{O}
  diffQ=q[i][0]-qaux[i][1]
   if diffQ != 0.0
      a[i][i]=(x[i][1]-olddx[i][1])/diffQ
  else
      a[i][i]=0.0
  end
  
  return nothing
end
#nupdateLinear differ from updateLinear in that nupdate uses with qminus instead of qaux=qminus+e*dq
#= function nupdateLinearApprox(i::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},qminus::Vector{Float64},olddx::Vector{MVector{O,Float64}})where{O}
  diffQ=q[i][0]-qminus[i]
    if diffQ != 0.0
     a[i][i]=(x[i][1]-olddx[i][1])/diffQ
  else
      a[i][i]=0.0
  end
  return nothing
end =#


function updateOtherApprox(k::Int,j::Int,x::Vector{Taylor0},q::Vector{Taylor0},a::Vector{Vector{Float64}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},simt)where{O}
  diffQ=q[j][0]-qaux[j][1]
  if diffQ != 0.0
    a[k][j]=(x[k][1]-olddx[k][1])/diffQ
  else
    a[k][j]=0.0
  end
  return nothing
end

