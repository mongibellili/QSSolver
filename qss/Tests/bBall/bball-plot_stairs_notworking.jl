

using formalqss
#= using Plots;
gr(); =#
function test()
    odeprob = @NLodeProblem begin
        parameter1=30.0# cache can be dynamic....parameters take this feature
        parameter2=1e+3
        u = [30.0,0.0]
        discrete = [25.0,0.0]#discrete1==stair; discrete2=contact
        du[1] =u[2]
        du[2] =-9.8-discrete[2]*((u[1]-discrete[1])*parameter2+parameter1*u[2])
       #=  du[3]=2.0*u[4]
        du[4]=-0.1*u[4] =#

        if discrete[1]-u[1]>0   #5*discrte gave error
            discrete[2]=1.0         
        else
            discrete[2]=0.0                                    
        end

        if u[2]>0
            discrete[1]=discrete[1]-5.0
       #=  else
            discrete[3]=5.0 =#
        end
        #= if u[3]-17.0+discrete[1]>0
            discrete[1]=discrete[1]-2.0
        else
            discrete[3]=5.0
        end =#
    end
    sol= QSS_Solve(odeprob,8.0,nmliqss2())
    save_Sol(sol)
end
test()
