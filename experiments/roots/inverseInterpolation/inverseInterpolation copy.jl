#@show cubic5(-1.33,2.0,1.0,6.0)#2.515556238254016 =#...cubic:18ns
#@btime allrealquadratic(1.33,1.0,-6.0)#(-2.5329305695911355, 1.7810508703430152)#1.647 ns (0 allocations: 0 bytes)
using BenchmarkTools
f(x)=1.33x*x+x-6.0
p0=-6.0
p1=f(1.0)#-3.67
p2=f(2.0)#1.3200000000000003
function quadratic(p0::F, p1::F, p2::F) where F <: AbstractFloat
    p123=p1+p2+p0
    psqr=p1*p1+p2*p2+p0*p0
    psqr2=p0*p0+p0*p0-p1-p2
    pdiff=p0+p0-p1-p2
    theta=-3.0(p0*p0-psqr/3.0)/psqr
    β=(theta+1.0)/((pdiff/psqr2)+p123/3.0-p0)
    α=(-3.0-β*pdiff)/psqr2
    γ=1.0-α*psqr/3.0-β*p123/3.0
end
#@btime quadratic(p0,p1,p2)#allocates
#@show quadratic(-6.0,-3.67,1.3200000000000003)#1.6687562560576625
#@btime quadratic(-6.0,-3.67,1.3200000000000003)
function quadratic(p0::F, p1::F, p2::F) where F <: AbstractFloat
    p123=p1+p2+p0
    psqr=p1*p1+p2*p2+p0*p0
    psqr2=p0*p0+p0*p0-p1-p2
    pdiff=p0+p0-p1-p2
    theta=-3.0(p0*p0-psqr/3.0)/psqr
    β=(theta+1.0)/((pdiff/psqr2)+p123/3.0-p0)
    α=(-3.0-β*pdiff)/psqr2
    γ=1.0-α*psqr/3.0-β*p123/3.0
end