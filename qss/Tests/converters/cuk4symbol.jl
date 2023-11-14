using qss
using BenchmarkTools


function test()
    odeprob = @NLodeProblem begin
      name=(cuk4,)
         C = 1e-4; L = 1e-4; R = 10.0;U = 24.0; T = 1e-4; DC = 0.25; ROn = 1e-5;ROff = 1e5;L1=1e-4;C1=1e-4;C2 = 1e-4;L2 = 1e-4;
         #discrete Rd(start=1e5), Rs(start=1e-5), nextT(start=T),lastT,diodeon;
         discrete = [1e5,1e-5,1e-4,0.0,0.0]
         Rd=discrete[1];Rs=discrete[2];nextT=discrete[3];lastT=discrete[4];diodeon=discrete[5]
         u[1:13]=0.0
         uc2=u[13]
         il2_1=u[i] ;il2_2=u[i-4] ;il2_3=u[i-8] ;il1_1=u[i+4] ;il1_2=u[i] ;il1_3=u[i-4] ;uc1_1=u[i+8];uc1_2=u[i+4] ;uc1_3=u[i] ;
         id1=(((il2_1+il1_1)*Rs-uc1_1)/(Rd+Rs))
         id2=(((il2_2+il1_2)*Rs-uc1_2)/(Rd+Rs))
         id3=(((il2_3+il1_3)*Rs-uc1_3)/(Rd+Rs))
        for i=1:4    #il2
          du[i] =(-uc2-Rs*id1)/L2
        end
        for i=5:8    #il1
          du[i]=(U-uc1_2-id2*Rs)/L1
        end
        for i=9:12    #uc1
          du[i]=(id3-il2_3)/C1
        end

        du[13]=(u[1]+u[2]+u[3]+u[4]-uc2/R)/C2



        if t-nextT>0.0 
          lastT=nextT
          nextT=nextT+T
          Rs=ROn
         
      end

      if t-lastT-DC*T>0.0 
          Rs=ROff
         
      end                          
      
  
      if diodeon*(((u[1]+u[5])*Rs-u[9])/(Rd+Rs))+(1.0-diodeon)*(((u[1]+u[5])*Rs-u[9])*Rd/(Rd+Rs))>0
        Rd=ROn
        diodeon=1.0
      else
        Rd=ROff
        diodeon=0.0
      end 

      if diodeon*(((u[2]+u[6])*Rs-u[10])/(Rd+Rs))+(1.0-diodeon)*(((u[2]+u[6])*Rs-u[10])*Rd/(Rd+Rs))>0
        Rd=ROn
        diodeon=1.0
      else
        Rd=ROff
        diodeon=0.0
      end 
           

      if diodeon*(((u[3]+u[7])*Rs-u[11])/(Rd+Rs))+(1.0-diodeon)*(((u[3]+u[7])*Rs-u[11])*Rd/(Rd+Rs))>0
        Rd=ROn
        diodeon=1.0
      else
        Rd=ROff
        diodeon=0.0
      end 

      if diodeon*(((u[4]+u[8])*Rs-u[12])/(Rd+Rs))+(1.0-diodeon)*(((u[4]+u[8])*Rs-u[12])*Rd/(Rd+Rs))>0
        Rd=ROn
        diodeon=1.0
      else
        Rd=ROff
        diodeon=0.0
      end 


    end
    @show odeprob
    tspan=(0.0,0.001)
   sol= solve(odeprob,nmliqss2(),abstol=1e-4,reltol=1e-3,tspan)
            
  save_Sol(sol)
  #save_Sol(sol,xlims=(1.9e-2,2e-2) ,ylims=(-6.0,8))
end
#@btime 
test()