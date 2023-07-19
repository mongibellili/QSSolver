using formalqss
#= using BSON
using BenchmarkTools
using XLSX =#
using TimerOutputs
using BenchmarkTools
function test()
   # BSON.@load "formalqss/ref_bson/solVectAdvection_N1000d01_Feagin14e-12.bson" solFeagin14VectorN1000d01
  prob=@NLodeProblem begin
    name=(adrN1000d01,)
    u[1:333]=1.0
    u[334:1000]=0.0
    _dx=100.0#1/dx=N/10=1000/10
    a=1.0
    d=0.1
    r=1000.0
   
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:999  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[1000]=-a*_dx*(u[1000]-u[999])+d*_dx*_dx*(2.0*u[999]-2.0*u[1000])+r*u[1000]*u[1000]*(1.0-u[1000]) 
    
end
    println("start solving")
   
    @btime solliqss2=QSS_Solve($prob,liqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    
end
test()
