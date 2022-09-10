#@show cubic5(-1.33,2.0,1.0,6.0)#2.515556238254016 =#...cubic:18ns
#@btime allrealquadratic(1.33,1.0,-6.0)#(-2.5329305695911355, 1.7810508703430152)#1.647 ns (0 allocations: 0 bytes)
using BenchmarkTools
#using LinearSolve
f(x)=1.33x*x+x-6.0
a=1.33
b=1.0
c=-6.0
#= p1=-6.0
p2=f(1.0)#-3.67
p3=f(2.0)#1.3200000000000003 =#

function quadratic(p0::Float64, p1::Float64, p2::Float64)::Float64
    p123=p1+p2+p0
    psqr=p1*p1+p2*p2+p0*p0
    psqr2=p0*p0+p0*p0-p1-p2
    pdiff=p0+p0-p1-p2
    theta=-3.0(p0*p0-psqr/3.0)/psqr
    β=(theta+1.0)/((pdiff/psqr2)+p123/3.0-p0)
    α=(-3.0-β*pdiff)/psqr2
    γ=1.0-α*psqr/3.0-β*p123/3.0
end
function lagrangeinversequadratic(x1::Float64,x2::Float64,x3::Float64,p1::Float64, p2::Float64, p3::Float64)::Float64
    p12=p1-p2
    p13=p1-p3
    p23=p2-p3
    finv=(p2*p3*x1*p23-p1*p3*x2*p13+p1*p2*x3*p12)/(p12*p13*p23)
end
function inversequadratic(x1::Float64,x2::Float64,x3::Float64,p1::Float64, p2::Float64, p3::Float64)::Float64
    sum=x1+x2+x3
    sum2=x1+x1-x2-x3
    psum=p1+p2+p3
    psqr=p1*p1+p2*p2+p3*p3
    psqr2=p1*p1+p1*p1-p2*p2-p3*p3
    psum2=p1+p1-p2-p3

    T=p3-psum2*p3*p3/psqr2+psqr*psum2/(3*psqr2)-psum/3
    R=sum2*p3*p3/psqr2+sum/3-psqr*sum2/(3*psqr2)
    β=(x3-R)/T
    α=(sum2-β*psum2)/psqr2
    γ=sum/3.0-α*psqr/3.0-β*psum/3.0
end
#@show quadratic(-6.0,-3.67,1.3200000000000003)#1.6687562560576625 when picked points are away
x1=1.5
x2=1.7
x3=1.95
p1=f(x1)
p2=f(x2)
p3=f(x3)
#= 
A=[p1*p1 p1 1;p2*p2 p2 1;p3*p3 p3 1;]

b=[x1,x2,x3] =#


function inversequadLinearSol(A::Matrix{Float64},b::Vector{Float64})
    prob=LinearProblem(A,b)
    sol=solve(prob)[3]
end
#@show inversequadLinearSol(A,b)
#@btime inversequadLinearSol(A,b)

#= C=[2.27256  -1.5075   1.0;
0.20821  -0.4563   1.0;
1.0147    1.00733  1.0]
vec=[1.5,1.7,1.95] =#
#@btime inversequadLinearSol(C,vec)

@btime lagrangeinversequadratic(x1,x2,x3,p1,p2,p3)#7.079 ns (0 allocations: 0 bytes)
#@show inversequadratic(x1,x2,x3,p1,p2,p3)#1.7814949781369704
#@btime inversequadratic(x1,x2,x3,p1,p2,p3)#35.285 ns (0 allocations: 0 bytes)
#@show inversequadratic(0.0,1.0,2.0,-6.0,-3.67,1.3200000000000003)#1.88
#@btime quadratic(p0,p1,p2)#27ns with 0 allocations
#@btime quadratic(-6.0,-3.67,1.3200000000000003)

function allreallinear(a::F, b::F) where F <: AbstractFloat
    if a == 0.0
      res = Inf
    else
      res = -b / a
    end
    return res
end
function allrealquadratic(a::F, b::F, c::F) where F <: AbstractFloat
    if a == 0.0
      res1=Inf
      res2 = allreallinear(b,c)
    else
      if c == 0.0
        res1=0.0
        res2 = allreallinear(a,b)
      else
        if b == 0.0
          r = -c / a
          if r <= 0.0
            res1=Inf
            res2=Inf
          else
            res1 = sqrt(r)
            res1 = -res1
          end
        else
          Δ = 1.0 - 4c*a / (b*b)
          if Δ < 0.0
            res1=Inf
            res2=Inf
          else
            q = -0.5*(1.0+sign(b)*sqrt(Δ))*b
            res1 = q / a
           
            res2=c / q
            
          end
        end
      end
    end
    return res1,res2
end
#@show allrealquadratic(1.33,1.0,-6.0)
#@btime allrealquadratic(a,b,c)#38.675 ns (1 allocation: 32 bytes)
#@btime allrealquadratic(1.33,1.0,-6.0)#1.907 ns (0 allocations: 0 bytes)