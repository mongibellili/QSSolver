using StaticArrays
using BenchmarkTools
#= macro make_simulator(settingstates)
    esc(quote struct ModelSettings
            finalTime:: Float64  
            initCond:: SVector{$settingstates,Float64}
            end
        struct ModelSimulator
            order::Int
            tx:: MVector{$settingstates,Float64}
            end
    end)
end
#**********user space1**********
states=3
initCond =  rand(states)
#display(initCond)
 @make_simulator (states )
#**********end user space1**********

function integrate(set::ModelSettings,sim::ModelSimulator)
#=     display(set);println()
    display(sim);println() =#
        for i = 1:1e+3 for j = 1:3
        sim.tx[j] =set.initCond[j] 
    end end
   
end
function simulate(initCond::Vector{Float64})
    states=3
    order=1
    tx = @MVector zeros(states)
    set=ModelSettings(2.0,initCond)
    sim=ModelSimulator(order,tx)
    integrate(set,sim)
    #display(tx)
end
#**********end user space2**********    
@btime simulate(initCond)
#**********end user space2********** =#



#--------------------------------------------------------------old approach----------------------------------------
struct ModelSettings{N}
            finalTime:: Float64  
            initCond:: SVector{N,Float64} 
            end

struct ModelSimulator{N}
            order::Int
            tx:: MVector{N,Float64}  
end
   

#**********user space1**********


#**********end user space1**********
function integrate(set::ModelSettings,sim::ModelSimulator)
#=     display(set);println()
    display(sim);println() =#
    for i = 1:1e+0 for j = 1:3
        sim.tx[j] =set.initCond[j] 
    end end
end
function simulate(initCond::SVector{N,Float64}) where {N}
    states=3
    order=1
    tx = @MVector zeros(states)
    set=ModelSettings(2.0,initCond)
    sim=ModelSimulator(order,tx)
    integrate(set,sim)
   # display(tx)
end
#**********end user space2**********  
states=3
initCond = @SVector rand(states)  
#display(initCond)
@btime simulate(initCond)
#**********end user space2**********