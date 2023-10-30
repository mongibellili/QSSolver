using qss
using BenchmarkTools
#= using Plots;
gr(); =#

function test()
    odeprob = @NLodeProblem begin
        #= parameter1=3000.0# cache can be dynamic....parameters take this feature
        parameter2=0.00001 =#
        L=1e3;R=10.0;uu=311.0;ROn = 1.0; ROff = 0.0;w=314.16
        u = [0.5,0.1]
       discrete = [1.0]
        du[1] =1000.0*(uu*sin(w*t)-u[1])
        du[2]=L*(u[1]-u[2]*(R))*discrete[1]


      


        if -u[1]>0.0 
            discrete[1]=ROff
            u[2]=0.0
        else
            discrete[1] = ROn
        end

                                     
        
       #=  if -u[1]>0.0 
            discrete[1] = ROn;
            
        end =#
           
    end
   sol= QSS_Solve(odeprob,qss2(),dQmin=1e-3,dQrel=1e-2,finalTime=0.06)
  # @show sol
 # @show 5
  save_Sol(sol)
  # save_Sol(sol,xlims=(0.0,0.06) ,ylims=(-2.04e-1,2.06e-1))
end
#@btime 
test()
