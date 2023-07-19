using formalqss

using TimerOutputs

function test()
 #reset_timer!(to)
  prob=@NLodeProblem begin
    name=(adrN10d01,)
    u[1:3]=1.0
    u[4:10]=0.0
    _dx=1.0#1/dx=N/10=10/10
    a=1.0
    d=0.1
    r=1000.0
   
    du[1] = -a*_dx*(u[1]-0.0)+d*_dx*_dx*(u[2]-2.0*u[1]+0.0)+r*u[1]*u[1]*(1.0-u[1]) 
    for k in 2:9  
        du[k]=-a*_dx*(u[k]-u[k-1])+d*_dx*_dx*(u[k+1]-2.0*u[k]+u[k-1])+r*u[k]*u[k]*(1.0-u[k]) ;
    end 
    du[10]=-a*_dx*(u[10]-u[9])+d*_dx*_dx*(2.0*u[9]-2.0*u[10])+r*u[10]*u[10]*(1.0-u[10]) 
    
end
    println("start solving")
   
 # solliqss2=QSS_Solve(prob,liqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
  #  solmliqss2=QSS_Solve(prob,mliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
  # solnliqss2=QSS_Solve(prob,nliqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
   solnmliqss2=QSS_Solve(prob,liqss2(),dQmin=1e-5,saveat=0.01,dQrel=1e-5,finalTime=10.0)#
    
  


  
  #pp=print_timer()
 # show(to)
#display(pp)
end
test()
