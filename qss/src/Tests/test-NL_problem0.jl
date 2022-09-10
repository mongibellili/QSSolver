using qss
using BenchmarkTools

using Plots;
gr();
function test()
    isPlot = true
    ssettings = SimSettings(5.0, saveat(0.1), qss2())###############################6.47235
    #ssettings = SimSettings()
    odeprob = @NLodeProblem begin
        cacheT=10
        u = [10.0,0.0]
        discrete = [0.0]

        du[1] =u[2]
        du[2] =-9.8-(discrete[1])*(1e+5*u[1]+30.0*u[2])
        #du[2] =-9.8-(discrete[1])*(1e+4 *(u[1]-discrete[2])+30.0*(u[2]+discrete[2]))*discrete[3]
       # dq[2] =-9.8-(d[1])*(300*q[1])
        if u[1]>0   #5*discrte gave error
            discrete[1]=0.0               
        else
            discrete[1]=1.0  
                                  
        end
#=         if 2.0-t>0
            discrete[2]=0.0
        else
            discrete[2]=5.0
            discrete[3]=0.001
           # q[2]=0.0 # still have not deal with elements of block
        end
        if 5.0-t>0
            #discrete[2]=0.0
        else
            discrete[2]=10.0
            discrete[3]=0.001
           # q[2]=0.0 # still have not deal with elements of block
        end =#
   
    end

#display(odeprob)
    sol = QSS_Solve(ssettings, odeprob)
    #display(sol)
    if isPlot



        temp1 = [] # in case the savedvars are of type taylor
        temp2 = []
       # temp3 = []
       # temp4 = []
        #if sol isa Tuple
            for i = 1:length(sol[2][1])
                push!(temp1, sol[2][1][i].coeffs[1])
                push!(temp2, sol[2][2][i].coeffs[1])
              #  push!(temp3, sol[2][1][i].coeffs[2]) #derx1 which is the same as x2
              #  push!(temp4, sol[2][2][i].coeffs[2])
            end

#= temp1=sol[2][1]
temp2=sol[2][2] =#

            labels = ["x1" "x2" "derx1" "derx2"]
            markercolors = [:green :orange :black :purple :red   :yellow :brown :white]
            display(plot!(sol[1], temp1,label="x1"))#,line=(:dot, 2),color=:green))
            display(plot!(sol[1], temp2,label="x2",line=(:dot, 4),color=:blue))
          #  display(plot!(sol[1], temp3,label="derx1",color=:orange))
           # display(plot!(sol[1], temp4,label = "derx2",color = :purple, xlims = (0.0,15.0),ylims = (-50.0,100.0)))
           # ylims!((-10,100))
           
            #, xlims = (0,0.00005),ylims = (1,1.0001)))
       # end
       println("done")
        readline()

    end




end

#@btime test()
test()

#= 
 using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
       #=  du[1] = u[2]#*cos(t) #*30*exp(t)
        du[2] = -u[1] - u[2]#+sqrt(t) # +1#+10-t*t+t +cos(t) =#
       # du[2]=(1 - u[1]^2) * u[2] - u[1] 
      # du[2]=1/(t+1) + sin(t)*sqrt(t)
     # du[2]=sqrt(t)


     du[1] = u[2]-u[1]+3.2
     #du[1] = 0.0+q[2]-0.0
     # inside I have T for both number and Taylor{T}...this requires to use 
                             #float cuz taylor used float....to be fixed later 
     du[2] =1.0-u[1]-u[2]


    
    end
    tspan = (0.0,5)
    u0 = [1.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
    println("done")
    readline()
end
#@btime odeDiffEquPackage()  #62.045 μs (439 allocations: 28.05 KiB) u2-u1+3.2  1-u1-u2
 odeDiffEquPackage() =#



































#=  using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
        du[1] = u[2]#*cos(t) #*30*exp(t)
        du[2] = -u[1] - u[2]#+sqrt(t) # +1#+10-t*t+t +cos(t)
       # du[2]=(1 - u[1]^2) * u[2] - u[1] 
      # du[2]=1/(t+1) + sin(t)*sqrt(t)
     # du[2]=sqrt(t)

    end
    tspan = (0.0,5)
    u0 = [1.0,2.4]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    display(plot!(sol,line=(:dot, 4)))
end
odeDiffEquPackage() =#



#@btime modqss.qss_Integrate(initCond,jacobian,order)


        #output=
#=         jacobian= StaticArrays.SVector{2, SymEngine.Basic}[[0, 1], [-100000.0*d1, -30.0*d1]]
        discJac= StaticArrays.SVector{1, SymEngine.Basic}[[0], [-100000.0*q1 - 30.0*q2]]
        zc_jac= StaticArrays.SVector{2, SymEngine.Basic}[[1, 0]]
        ZC_discJac= StaticArrays.SVector{1, SymEngine.Basic}[[0]]
        evHandlr= qss.EventHandlerStruct[qss.EventHandlerStruct{2, 1}(1, [NaN, NaN], [0.0]), qss.EventHandlerStruct{2, 1}(2, [NaN, NaN], [1.0])]
        SD= StaticArrays.SVector{2, Int64}[[0, 2],[1, 2] ]
        SZ= StaticArrays.SVector{1, Int64}[[1], [0]]
        HZ= StaticArrays.SVector{1, Int64}[[0], [0]]
        HD= StaticArrays.SVector{2, Int64}[[0, 2], [0, 2]] =#