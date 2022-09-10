
using BenchmarkTools
using qss
#= using Plots;
gr(); =#
function test()
    isPlot = false
    ssettings = SimSettings(5.0, saveat(0.1), qss2())###############################
    #ssettings = SimSettings()
    odeprob = @NLodeProblem begin
        cacheT=10 #recommend the cache to match the largest equation you have using the specs of "caches/ops"
        u = [1.0, 0.0]
        discrete = [0.0]
        du[1] = u[2]
        #du[1] = 0.0+q[2]-0.0
        # inside I have T for both number and Taylor{T}...this requires to use 
                                #float cuz taylor used float....to be fixed later 
        du[2] =-u[1]-u[2]

       #=  if p[1] >0   #still have not added 'usercode checking'....throw error msg
            d[1]=0.0
            #u[2]=3.3 #this requires the function to be altered inside the macro u[2].coeffs[1]=3.3
        else
            d[1]=1.0          
        end =#
        #= du[1] = u[2]-u[1]+3.2
        du[2] =1.0-u[1]-u[2] =#
       


    end

    #display(odeprob)
    sol = QSS_Solve(ssettings, odeprob)
    #display(sol[3]) # 320 while loops (steps) for ft=5 qss2 for 2 states...recompute called about 1.5times for this prob each step...
  #=   if isPlot
        
        #= temp1 = []
        temp2 = []
        #if sol isa Tuple
            for i = 1:length(sol[2][1])
                push!(temp1, sol[2][1][i].coeffs[1])
                push!(temp2, sol[2][2][i].coeffs[1])
            end
            display(plot!(sol[1], temp1))
            display(plot!(sol[1], temp2))
       # end =#
     
    display(plot!(sol[1], sol[2][1]))
    display(plot!(sol[1], sol[2][2]))

       println("done")
        readline()

    end =#




end

@btime test()
#display(@benchmark test())
#test()

#= 
 using OrdinaryDiffEq
function odeDiffEquPackage()
    function funcName(du,u,p,t)# api requires four args
       #=  du[1] = u[2]#*cos(t) #*30*exp(t)
        du[2] = -u[1] - u[2]#+sqrt(t) # +1#+10-t*t+t +cos(t) =#
       # du[2]=(1 - u[1]^2) * u[2] - u[1] 
      # du[2]=1/(t+1) + sin(t)*sqrt(t)
     # du[2]=sqrt(t)


     du[1] = u[2]
     #du[1] = 0.0+q[2]-0.0
     # inside I have T for both number and Taylor{T}...this requires to use 
                             #float cuz taylor used float....to be fixed later 
     du[2] =-u[1]-u[2]


    
    end
    tspan = (0.0,5)
    u0 = [1.0,0.0]
    prob = ODEProblem(funcName,u0,tspan)
    sol = solve(prob,BS3(),abstol = 1e-6, reltol = 1e-3)
   # display(sol)
    #display(plot!(sol,line=(:dot, 4)))
   # println("done")
   # readline()
end
#display(@benchmark odeDiffEquPackage()) #62.045 μs (439 allocations: 28.05 KiB) u2-u1+3.2  1-u1-u2
 #odeDiffEquPackage()
 @btime odeDiffEquPackage() =#



































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