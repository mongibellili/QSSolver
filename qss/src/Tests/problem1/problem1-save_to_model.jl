using qss

function test()
    odeprob = @NLodeProblem begin
       #=  u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        du[2] =-u[1]-u[2] =#



        #= u = [-4.0, 4.0]
        discrete = [0.0]
        du[1] = -u[1]-10.0*u[2]+0.2
        du[2] =10.0*u[1]-u[2]+1.2 =#
#------------twoVarSys15----------
       #=  u = [0.5, 0.0]
        discrete = [0.0]
        du[1] = -4.0*u[1]+0.5*u[2]
        du[2] =u[1]-0.25*u[2]  =#
        u = [1.0, 2.3]
        discrete = [0.0]
        du[1] = -1.01*u[1]-100.0*u[2]+1000.2
        du[2] =10.1*u[1]-1.001*u[2]+0.2



    end
  
    save_prob_to_model(odeprob,"/home/unknown/QS_Solver/qss/src/models/classicProblem.jl","twoVarSys16") #any location you want
   
end
test()
