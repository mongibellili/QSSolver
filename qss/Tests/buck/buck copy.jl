using qss
using BenchmarkTools
#= using Plots;
gr(); =#

function test()
    odeprob = @NLodeProblem begin
       
         C = 1e-4; L = 1e-4; R = 10;U = 24.0; T = 1e-4; DC = 0.5; ROn = 1e-5;ROff = 1e5;
         #discrete Rd(start=1e5), Rs(start=1e-5), nextT(start=T),lastT,diodeon;
         
       
         u = [0.0,0.0]
      
        du[1] =(-((u[1]-24.0))-u[2])/L
        du[2]=(u[1]-u[2]/R)/C

      
     #=  if -(discrete[5]*((u[1]*discrete[2]-U)/(discrete[1]+discrete[2]))+(1-discrete[5])*((u[1]*discrete[2]-U)*discrete[1]/(discrete[1]+discrete[2])))>0
        #= discrete[1]=ROn
        discrete[5]=1.0
      else =#
        discrete[1]=ROff
        discrete[5]=0.0
      end =#
           
    end
   sol= QSS_Solve(odeprob,nmliqss2(),dQmin=1e-4,dQrel=1e-3,finalTime=10.0)
   #@show sol
 # @show 5
  save_Sol(sol)
  # save_Sol(sol,xlims=(0.0,0.06) ,ylims=(-2.04e-1,2.06e-1))
end
#@btime 
test()
